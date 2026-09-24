import os
from pathlib import Path
from pathlib import Path
import numpy as np
import pymatgen.core as pmg
from pymatgen.analysis.graphs import StructureGraph, MoleculeGraph
from pymatgen.analysis.local_env import JmolNN
import networkx as nx

# Heteroatoms carried across a severed linkage so that BOTH building blocks
# keep the functional group at their connection point. Carbon and hydrogen are
# excluded: carrying those would grow the blocks into each other's backbone
# rather than preserving peripheral functionality.
_CARRYOVER_ELEMENTS = {"O", "N", "S", "P", "B", "F", "Cl", "Br", "I"}

def _count_linkages(cut_bonds, linkage_of=None):
    """Count distinct chemical linkages from a component's severed bonds.

    `cut_bonds` are (attachment_atom, external_partner) pairs. Two severed
    bonds belong to the same linkage when they share an endpoint, or when the
    cleavage rule that produced them tagged them with the same linkage id in
    `linkage_of` (keyed by frozenset({u, v})).

    The explicit tag is required for ring linkages whose severed bonds share
    no endpoint at all: a 1,4-dioxin is cut at `C2-O1` and `C3-O4`, which have
    four distinct endpoints yet form ONE connection. Endpoint sharing alone
    covers the benzoxazole case (both bonds meet at the C2 carbon) and leaves
    single-bond linkages (imine, boroxine) counted one apiece.
    """
    if not cut_bonds:
        return 0
    lg = nx.Graph()
    lg.add_nodes_from(range(len(cut_bonds)))
    for a in range(len(cut_bonds)):
        for b in range(a + 1, len(cut_bonds)):
            same_endpoint = (cut_bonds[a][0] == cut_bonds[b][0]
                             or cut_bonds[a][1] == cut_bonds[b][1])
            same_tag = False
            if linkage_of:
                ta = linkage_of.get(frozenset(cut_bonds[a]))
                tb = linkage_of.get(frozenset(cut_bonds[b]))
                same_tag = ta is not None and ta == tb
            if same_endpoint or same_tag:
                lg.add_edge(a, b)
    return nx.number_connected_components(lg)


class SBU:
    """Represents a Secondary Building Unit (Node or Linker) in a COF."""
    def __init__(self, molecule, indices=None):
        self.molecule = molecule # pymatgen.core.Molecule
        self.indices = indices
        self.smiles = "" # Placeholder for SMILES, in practice calculated via RDKit/OpenBabel
        self.molar_mass = molecule.composition.weight if molecule else 0

class FragmentationResult:
    """Holds the resulting nodes and linkers from COF fragmentation.

    `n_cut_bonds` and `linkage_types` record what the cleavage rules actually
    recognised. Zero severed bonds means NO known linkage chemistry matched, in
    which case the "blocks" below are just the structure's connected components
    (often still-periodic networks) rather than real building blocks. Callers
    report that instead of failing silently - it is the signature of a COF whose
    chemistry is missing from the rule set.
    """
    def __init__(self, nodes, linkers, n_cut_bonds=None, linkage_types=None):
        self.nodes = nodes
        self.linkers = linkers
        self.n_cut_bonds = n_cut_bonds
        self.linkage_types = linkage_types or []

class COF:
    """
    Main class for Covalent Organic Framework representation.
    Inspired by moffragmentor's MOF class.
    """
    def __init__(self, structure):
        self.structure = structure
        
    @classmethod
    def from_cif(cls, cif_path):
        """Initialize a COF object from a CIF file."""
        structure = pmg.Structure.from_file(str(cif_path))
        return cls(structure)
        
    def fragment(self) -> FragmentationResult:
        """
        Fragments the COF into nodes and linkers.
        
        In COFs, this requires cleaving specific dynamic covalent bonds 
        (e.g., Imine C=N, Boroxine B-O, Hydrazone, Azine, etc.).
        """
        # 1. Build a structure graph to represent the connectivity of the COF
        nn_strategy = JmolNN()
        sg = StructureGraph.with_local_env_strategy(self.structure, nn_strategy)
        
        # 2. Identify the bonds to cleave (The "Linkage" bonds)
        # ---------------------------------------------------------
        # THIS IS WHERE YOUR SPECIFIC COF CLEAVAGE LOGIC GOES
        # For a practical COF fragmenter, you would search for specific 
        # atomic motifs. For example, finding Carbon-Nitrogen double bonds
        # and removing those edges from the graph.
        # ---------------------------------------------------------
        edges_to_remove = []
        # Bonds severed by one multi-bond linkage share a tag here, keyed by
        # frozenset({u, v}); see _count_linkages.
        linkage_of = {}
        # Create an undirected version for component analysis
        undirected_graph = sg.graph.to_undirected()
        
        # Hydrogen is monovalent, but JmolNN scores by distance alone and will
        # perceive a short non-covalent contact as a bond. 526 is the case in
        # point: each of its six beta-ketoenamine hydrogens sits 1.09 A from its
        # own carbon AND 1.28 A from a neighbouring keto oxygen - the resonance
        # assisted H-bond drawn as a dashed line in any picture of a TpPa COF -
        # and the graph took both. Those six phantom edges bridge node and
        # linker, so cutting all six real linkages still left ONE 114-atom
        # component: no blocks at all, and the structure fell back to Path A
        # radius truncation. Dropping them recovers exactly the expected
        # decomposition, 3 x C12H8N4O4 linkers and 2 x C9H3O3 nodes.
        #
        # Keep only each hydrogen's nearest neighbour, which is always its
        # covalent parent, and drop the rest.
        for _h in [
            i for i in range(len(self.structure))
            if self.structure[i].specie.symbol == 'H'
        ]:
            _nbrs = list(undirected_graph.neighbors(_h))
            if len(_nbrs) <= 1:
                continue
            # Prefer the nearest HEAVY neighbour, not the nearest neighbour of
            # any kind. 70 has a disordered hydrogen cluster where H[23] sits
            # 1.07 A from its carbon but only 0.86 A from another H; picking
            # the bare minimum would keep that H-H pair and throw away the real
            # C-H bond.
            _heavy = [
                j for j in _nbrs
                if self.structure[j].specie.symbol != 'H'
            ]
            _pool = _heavy or _nbrs
            _keep = min(_pool, key=lambda j: self.structure.get_distance(_h, j))
            for _nb in _nbrs:
                if _nb == _keep:
                    continue
                while undirected_graph.has_edge(_h, _nb):
                    undirected_graph.remove_edge(_h, _nb)

        # Create a copy to preserve original connections for neighbor analysis
        original_undirected = undirected_graph.copy()
        
        # Helper: count heavy (non-H) neighbors of a node
        def heavy_degree(node):
            return sum(1 for nb in undirected_graph.neighbors(node)
                       if self.structure[nb].specie.symbol != 'H')

        # Helper: count H neighbours of a node
        def h_count(node):
            return sum(1 for nb in undirected_graph.neighbors(node)
                       if self.structure[nb].specie.symbol == 'H')

        # Helper: distance between two bonded sites
        def bond_length(u, v):
            return float(self.structure.get_distance(u, v))

        # Alkene carbons paired across a severed vinylene linkage. Carbon is
        # normally excluded from carry-over, but the partner of a cut C=C is
        # kept so each block retains an intact, planar sp2 vinyl terminus
        # (Ar-CH=CH2) instead of collapsing to a tetrahedral Ar-CH3.
        vinylene_partner = {}

        # Helper: check if an edge is in a small ring (e.g. porphyrin pyrrole ring)
        def in_small_ring(u, v, max_len=6):
            G_temp = undirected_graph.copy()
            G_temp.remove_edge(u, v)
            try:
                path = nx.shortest_path(G_temp, u, v)
                return len(path) <= max_len
            except nx.NetworkXNoPath:
                return False

        # Pyrazine linkages: two aromatic units fused through a pyrazine ring
        # (quinoxaline, phenazine, hexaazatriphenylene). The ring is the
        # linkage - it is what the condensation of an o-diamine with an
        # o-diketone forms - but every one of its bonds sits inside a ring, so
        # no rule here touched it and 662 and 663 came out as a patch of sheet
        # with no node or linker at all.
        #
        # Both nitrogens are severed from both of their carbons, which leaves
        # each nitrogen alone and the two carbon units separate. The carry-over
        # below then gives each unit its own copy of the nitrogens it was
        # bonded to, so 662 decomposes into the hexa-substituted benzene node
        # carrying six N and the tetra-substituted linker carrying four - the
        # two halves of the condensation, each showing its linkage
        # environment, exactly as for every other linkage type here.
        #
        # The test: a six-ring holding exactly two nitrogens para to each
        # other, each with exactly two heavy neighbours, both carbons, both in
        # that ring. Pyridine-type ring N (one per ring) and fused imidazoles
        # do not match, and neither does a pyrazine whose nitrogen carries a
        # substituent.
        pyrazine_bonds = {}
        pyrazine_bridge_atoms = set()

        def _heavy_nbrs_of(idx):
            return [nb for nb in undirected_graph.neighbors(idx)
                    if self.structure[nb].specie.symbol != 'H']

        def _is_bridge_n(idx):
            if self.structure[idx].specie.symbol != 'N':
                return None
            heavy = _heavy_nbrs_of(idx)
            if len(heavy) != 2:
                return None
            if any(self.structure[nb].specie.symbol != 'C' for nb in heavy):
                return None
            return heavy

        # Walk the ring by hand rather than through a cycle basis: a basis is
        # free to return one big cycle instead of the six-rings it is made of,
        # and the pyrazines of 662 come back inside a sixteen-membered cycle.
        for _n1 in list(undirected_graph.nodes()):
            _pair1 = _is_bridge_n(_n1)
            if _pair1 is None:
                continue
            _cA, _cB = _pair1
            found = None
            for _cA2 in _heavy_nbrs_of(_cA):
                if _cA2 == _n1 or self.structure[_cA2].specie.symbol != 'C':
                    continue
                for _n2 in _heavy_nbrs_of(_cA2):
                    if _n2 == _cA or _n2 == _n1:
                        continue
                    _pair2 = _is_bridge_n(_n2)
                    if _pair2 is None:
                        continue
                    _other = [c for c in _pair2 if c != _cA2]
                    if len(_other) != 1:
                        continue
                    _cB2 = _other[0]
                    if _cB2 == _cB or not undirected_graph.has_edge(_cB2, _cB):
                        continue
                    found = (_n2, _cA2, _cB2)
                    break
                if found:
                    break
            if not found:
                continue
            _n2 = found[0]
            # All four bonds of one pyrazine carry the same linkage id, keyed
            # on its two nitrogens. Without that they count as four separate
            # linkages, and a linker sitting between two pyrazines scores four
            # - over the three that mark a node - so 662 came back as five
            # nodes and no linker at all.
            _tag = ('pyrazine', frozenset((_n1, _n2)))
            for _n in (_n1, _n2):
                pyrazine_bridge_atoms.add(_n)
                for nb in _heavy_nbrs_of(_n):
                    if self.structure[nb].specie.symbol == 'C':
                        pyrazine_bonds[frozenset((_n, nb))] = _tag

        # Nitrogen-nitrogen linkages: azine (Ar-CH=N-N=CH-Ar), acylhydrazone
        # (Ar-CH=N-NH-CO-Ar), azo and hydrazo (Ar-N=N-Ar, Ar-NH-NH-Ar). In all
        # of them the N-N bond IS the linkage, and it is the bond to sever.
        #
        # Cutting the flanking C=N instead - which the imine rule below does
        # wherever that carbon has two heavy neighbours - leaves the whole
        # N-N unit dangling off one block, so that block terminates as =N-NH2
        # while the other side ends as a bare aryl. 1223 was the case in
        # point: all six of its bridges were cut at the C=N, and the fragment
        # came back with two nitrogens at some edges and one at others.
        # Severing the N-N puts one nitrogen on each side, and every edge
        # terminates the same way.
        #
        # The predicate is the bridging geometry, not any one chemistry: both
        # nitrogens must have exactly two heavy neighbours, one carbon and one
        # nitrogen, and the bond must lie outside any small ring. That covers
        # every N-N bridge in the CoRE-COF set - 34 azine, 66 acylhydrazone,
        # 15 azo/hydrazo, 2 mixed - while excluding ring systems (pyrazole,
        # pyridazine, tetrazine) and pendant hydrazides, whose terminal
        # nitrogen has no carbon of its own.
        #
        # The carbon on either side keeps its own chemistry, so the caps come
        # out right per family without a special case: an azine gives two
        # Ar-CH=NH aldimines, an acylhydrazone gives one aldimine and one
        # primary amide, an azo or hydrazo bridge gives two anilines.
        azine_nn_bonds = set()
        azine_nitrogens = set()
        for _u, _v in undirected_graph.edges():
            if (self.structure[_u].specie.symbol != 'N'
                    or self.structure[_v].specie.symbol != 'N'):
                continue
            _ok = True
            for _n in (_u, _v):
                _heavy = [nb for nb in undirected_graph.neighbors(_n)
                          if self.structure[nb].specie.symbol != 'H']
                if len(_heavy) != 2:
                    _ok = False
                    break
                if sorted(self.structure[nb].specie.symbol
                          for nb in _heavy) != ['C', 'N']:
                    _ok = False
                    break
            if _ok and not in_small_ring(_u, _v):
                azine_nn_bonds.add(frozenset((_u, _v)))
                azine_nitrogens.update((_u, _v))

        for u, v, data in undirected_graph.edges(data=True):
            atom_u = self.structure[u].specie.symbol
            atom_v = self.structure[v].specie.symbol
            
            bond_pair = set([atom_u, atom_v])
            
            # Heuristic: break common COF forming bonds.
            if bond_pair == {'B', 'O'}:
                edges_to_remove.append((u, v))
                linkage_of[frozenset((u, v))] = ('boroxine', frozenset((u, v)))
            elif bond_pair == {'C', 'N'} and frozenset((u, v)) in pyrazine_bonds:
                edges_to_remove.append((u, v))
                linkage_of[frozenset((u, v))] = pyrazine_bonds[frozenset((u, v))]
            elif bond_pair == {'C', 'N'}:
                # Cut the imine C=N double bond itself (where Carbon has heavy degree == 2)
                # so that the terminal Nitrogen atoms stay with the amine-derived building block
                # and terminal Oxygen/Carbon atoms stay with the aldehyde-derived building block.
                # Protect C-N bonds that are part of small rings (porphyrin, triazine, pyridine).
                c_idx = u if atom_u == 'C' else v
                n_idx = v if c_idx == u else u
                # A linkage nitrogen always bridges two units, so it has two
                # heavy neighbours. A TERMINAL nitrogen is a substituent - a
                # nitrile (Ar-C#N) or a free amine - and cutting it severs a
                # functional group instead of a linkage. 760.cif is the case in
                # point: its eight nitriles were being cut, which decomposed
                # nothing (the framework stayed one 236-atom component) but was
                # enough to make `edges_to_remove` non-empty, so the biaryl
                # fallback never ran and the structure dropped to Path D.
                if (
                    heavy_degree(c_idx) == 2
                    and heavy_degree(n_idx) >= 2
                    and n_idx not in azine_nitrogens
                    and not in_small_ring(u, v)
                ):
                    edges_to_remove.append((u, v))
                    linkage_of[frozenset((u, v))] = ('imine', frozenset((u, v)))
            elif bond_pair == {'C'}:
                # Vinylene (sp2-carbon) COF linkage: Ar-CH=CH-Ar, formed by
                # Knoevenagel/aldol condensation. Both alkene carbons carry
                # exactly one H and exactly two heavy neighbours (each other
                # plus one aryl carbon), the C=C is short (~1.34 A) while the
                # flanking Ar-C bonds are long (~1.48 A), and the bond lies
                # outside any small ring - which is what separates it from an
                # ordinary aromatic CH=CH edge, whose atoms look identical
                # locally. Without this rule a vinylene COF has no recognised
                # linkage at all: nothing is cut, the "blocks" come back as
                # whole periodic networks, and the caller drops to its crudest
                # fallback path.
                if (
                    heavy_degree(u) == 2
                    and heavy_degree(v) == 2
                    and h_count(u) == 1
                    and h_count(v) == 1
                    and bond_length(u, v) <= 1.42
                    and not in_small_ring(u, v)
                ):
                    edges_to_remove.append((u, v))
                    linkage_of[frozenset((u, v))] = ('vinylene', frozenset((u, v)))
                    vinylene_partner[u] = v
                    vinylene_partner[v] = u
            elif bond_pair == {'N'}:
                if frozenset((u, v)) in azine_nn_bonds:
                    edges_to_remove.append((u, v))
                    linkage_of[frozenset((u, v))] = ('N-N bridge', frozenset((u, v)))

        # Benzoxazole / oxazole linkage. The C2 carbon of the five-membered
        # oxazole ring is bonded to BOTH the ring oxygen and the ring nitrogen
        # and is the former aldehyde carbon; severing both of those bonds
        # separates the aldehyde-derived building block from the
        # aminophenol-derived one. The C-N bond here sits inside a small ring,
        # so the imine rule above deliberately skips it and it must be handled
        # explicitly.
        for c_idx in list(undirected_graph.nodes()):
            if self.structure[c_idx].specie.symbol != 'C':
                continue
            o_nbs = [nb for nb in undirected_graph.neighbors(c_idx)
                     if self.structure[nb].specie.symbol == 'O']
            n_nbs = [nb for nb in undirected_graph.neighbors(c_idx)
                     if self.structure[nb].specie.symbol == 'N']
            if len(o_nbs) != 1 or len(n_nbs) != 1:
                continue
            o_idx, n_idx = o_nbs[0], n_nbs[0]
            probe = undirected_graph.copy()
            probe.remove_node(c_idx)
            try:
                # O-x-y-N (4 nodes) closes a five-membered ring through c_idx.
                if len(nx.shortest_path(probe, o_idx, n_idx)) == 4:
                    edges_to_remove.append((c_idx, o_idx))
                    edges_to_remove.append((c_idx, n_idx))
                    tag = ('oxazole', c_idx)
                    linkage_of[frozenset((c_idx, o_idx))] = tag
                    linkage_of[frozenset((c_idx, n_idx))] = tag
            except nx.NetworkXNoPath:
                pass

        # 1,4-dioxin linkage (dioxin-linked COFs, e.g. HHTP + perfluoroarene).
        # The six-membered ring O1-C2-C3-O4-C5-C6 fuses two different aryl
        # systems: {C2,C3} belong to one, {C5,C6} to the other. Severing the
        # two C-O bonds on ONE side separates them. The oxygens are kept with
        # the LARGER aryl system, which is the polyol-derived building block
        # (hexahydroxytriphenylene in the canonical case).
        o_atoms = [i for i in undirected_graph.nodes()
                   if self.structure[i].specie.symbol == 'O']
        if o_atoms:
            simple = nx.Graph(undirected_graph)
            no_o = simple.copy()
            no_o.remove_nodes_from(o_atoms)
            aryl_block = {}
            for comp in nx.connected_components(no_o):
                rep, size = min(comp), len(comp)
                for a in comp:
                    aryl_block[a] = (rep, size)

            seen_rings = set()
            for o1 in o_atoms:
                c_nbs = [nb for nb in simple.neighbors(o1)
                         if self.structure[nb].specie.symbol == 'C']
                if len(c_nbs) != 2:
                    continue
                for a in c_nbs:
                    for b in simple.neighbors(a):
                        if b == o1 or self.structure[b].specie.symbol != 'C':
                            continue
                        for o4 in simple.neighbors(b):
                            if o4 == a or self.structure[o4].specie.symbol != 'O':
                                continue
                            # O1-a-b-O4 found; require the ring to close back
                            # through a second, disjoint C-C bridge.
                            other_a = [x for x in c_nbs if x != a]
                            o4_cs = [x for x in simple.neighbors(o4)
                                     if self.structure[x].specie.symbol == 'C' and x != b]
                            if not other_a or not o4_cs:
                                continue
                            c6 = other_a[0]
                            if not any(simple.has_edge(c6, c5) for c5 in o4_cs):
                                continue
                            ring_key = frozenset((o1, o4, a, b, c6))
                            if ring_key in seen_rings:
                                continue
                            if a not in aryl_block or b not in aryl_block:
                                continue
                            if aryl_block[a][0] != aryl_block[b][0]:
                                continue  # a and b must be the same aryl side
                            c5 = next(c5 for c5 in o4_cs if simple.has_edge(c6, c5))
                            if c5 not in aryl_block or c6 not in aryl_block:
                                continue
                            if aryl_block[c5][0] != aryl_block[c6][0]:
                                continue
                            if aryl_block[a][0] == aryl_block[c5][0]:
                                continue  # both sides in one system: not a linkage
                            seen_rings.add(ring_key)
                            # cut the O-C bonds on the smaller aryl side
                            if aryl_block[a][1] <= aryl_block[c5][1]:
                                pair = ((o1, a), (o4, b))
                            else:
                                pair = ((o1, c6), (o4, c5))
                            tag = ('dioxin', tuple(sorted(ring_key)))
                            for u, v in pair:
                                edges_to_remove.append((u, v))
                                linkage_of[frozenset((u, v))] = tag

        # Biaryl fallback for all-hydrocarbon COFs. Frameworks such as 463.cif
        # (C48H30) are built purely from aromatic rings joined by direct
        # aryl-aryl C-C single bonds, so none of the linkage chemistries above
        # matches and NOTHING is cut: `fragment()` then returns whole periodic
        # networks (or nothing at all) and the caller drops to its crudest
        # fallback path. This rule runs ONLY when no other linkage was found,
        # so any structure with recognised chemistry is completely unaffected -
        # important, because most COFs contain biaryl bonds inside their
        # linkers that must NOT be cut.
        #
        # Which biaryl bonds to sever is decided by topology, not by length.
        # Rings are fused into ring systems and the systems are connected by
        # the biaryl bonds; a system reached by >=3 such bonds is a net vertex
        # (a node), while 2-connected systems are pieces of a strut. Severing
        # only the bonds that touch a node therefore keeps a two-ring strut
        # intact as one biphenyl linker instead of splitting it into two bare
        # rings. On 463 this cuts exactly the six node-linker bonds and keeps
        # the three intra-biphenyl ones, purely from connectivity.
        if not edges_to_remove:
            simple = nx.Graph(undirected_graph)
            small_rings = [set(r) for r in nx.minimum_cycle_basis(simple) if len(r) <= 7]
            if small_rings:
                fuse = nx.Graph()
                fuse.add_nodes_from(range(len(small_rings)))
                for a in range(len(small_rings)):
                    for b in range(a + 1, len(small_rings)):
                        if small_rings[a] & small_rings[b]:
                            fuse.add_edge(a, b)
                system_of = {}
                for k, comp in enumerate(nx.connected_components(fuse)):
                    for r in comp:
                        for atom in small_rings[r]:
                            system_of[atom] = k

                biaryl = []
                for u, v in simple.edges():
                    if (
                        self.structure[u].specie.symbol == 'C'
                        and self.structure[v].specie.symbol == 'C'
                        and u in system_of
                        and v in system_of
                        and system_of[u] != system_of[v]
                        and not in_small_ring(u, v)
                    ):
                        biaryl.append((u, v))

                degree = {}
                for u, v in biaryl:
                    degree[system_of[u]] = degree.get(system_of[u], 0) + 1
                    degree[system_of[v]] = degree.get(system_of[v], 0) + 1
                # Without a >=3-connected ring system there is no net vertex to
                # anchor the deconstruction, so leave the structure alone.
                if any(d >= 3 for d in degree.values()):
                    for u, v in biaryl:
                        if degree.get(system_of[u], 0) >= 3 or degree.get(system_of[v], 0) >= 3:
                            edges_to_remove.append((u, v))
                            linkage_of[frozenset((u, v))] = ('biaryl', frozenset((u, v)))

        # 2b. Orphan guard. A linkage rule matches a bond pattern, not a whole
        # linkage, so it can sever EVERY bond around a small bridging unit and
        # leave it floating: an Ar-NH-Ar secondary amine loses both C-N bonds
        # (931), an -N=N-N- triazene chain is cut on all sides (645), a C
        # bridging two N likewise (1226). The residue is a 1-3 heavy-atom
        # "linker" that is meaningless as a QM fragment. Restore one cut per
        # orphan so the unit stays attached to its largest neighbour and becomes
        # a proper terminus instead (Ar-NH- then caps to Ar-NH2).
        MIN_STRUT_HEAVY = 4
        for _ in range(12):
            trial = nx.Graph(undirected_graph)
            trial.remove_edges_from(edges_to_remove)
            comps = list(nx.connected_components(trial))
            if len(comps) <= 1:
                break
            comp_of, heavy_of = {}, {}
            for k, comp in enumerate(comps):
                heavy_of[k] = sum(
                    1 for i in comp if self.structure[i].specie.symbol != 'H'
                )
                for i in comp:
                    comp_of[i] = k
            # A pyrazine nitrogen is meant to be left on its own: it is
            # carried into both blocks below, not exported as a block, so
            # restoring one of its bonds here would undo the linkage cut.
            orphans = [
                k for k, h in heavy_of.items()
                if h < MIN_STRUT_HEAVY
                and not all(
                    i in pyrazine_bridge_atoms
                    or self.structure[i].specie.symbol == 'H'
                    for i in comps[k]
                )
            ]
            if not orphans:
                break
            restored = False
            for k in orphans:
                # Cuts with exactly one endpoint inside this orphan.
                cands = [
                    (u, v) for (u, v) in edges_to_remove
                    if (comp_of.get(u) == k) != (comp_of.get(v) == k)
                ]
                if not cands:
                    # Nothing was severed here: it is a guest/solvent molecule
                    # that was already disconnected in the parent, not an
                    # orphan we created. Leave it for the zero-cut skip below.
                    continue
                best = max(
                    cands,
                    key=lambda e: heavy_of.get(
                        comp_of.get(e[1]) if comp_of.get(e[0]) == k else comp_of.get(e[0]), 0
                    ),
                )
                edges_to_remove.remove(best)
                linkage_of.pop(frozenset(best), None)
                restored = True
            if not restored:
                break

        # 3. Cleave the bonds in the graph
        undirected_graph.remove_edges_from(edges_to_remove)
        
        # Helper: Create a set of cut edges for fast lookup
        cut_edges = set()
        for u, v in edges_to_remove:
            cut_edges.add((u, v))
            cut_edges.add((v, u))
        
        # 4. Extract disconnected subgraphs (these are your fragments)
        nodes = []
        linkers = []
        
        # Helper: Map atom to component for topological coordination calculation
        components = list(nx.connected_components(undirected_graph))
        atom_to_comp = {}
        for comp_id, comp in enumerate(components):
            for atom in comp:
                atom_to_comp[atom] = comp_id
                
        # Iterating over connected components
        for comp_id, comp in enumerate(components):
            comp_indices = list(comp)
            
            # Skip if it's the whole unfragmented structure or a single atom
            if len(comp_indices) == len(self.structure) or len(comp_indices) <= 1:
                continue
            
            # Determine connection points based on the number of attachment atoms
            # (atoms in this component that have bonds that were cut)
            attachment_atoms = set()
            cut_bonds = []

            # BFS to unwrap coordinates across periodic boundaries
            unwrapped_coords = {}
            start_node = comp_indices[0]
            unwrapped_coords[start_node] = self.structure[start_node].coords
            
            queue = [start_node]
            visited = {start_node}
            
            while queue:
                u = queue.pop(0)
                u_frac_unwrapped = self.structure.lattice.get_fractional_coords(unwrapped_coords[u])
                for v in undirected_graph.neighbors(u):
                    if v in comp_indices and v not in visited:
                        _, image = self.structure.lattice.get_distance_and_image(
                            u_frac_unwrapped, self.structure[v].frac_coords
                        )
                        v_frac = self.structure[v].frac_coords + image
                        unwrapped_coords[v] = self.structure.lattice.get_cartesian_coords(v_frac)
                        visited.add(v)
                        queue.append(v)
            
            sites = []
            site_indices = []
            carried = []

            for i in comp_indices:
                # Use the unwrapped coordinates instead of raw unit cell coordinates
                site_coord = unwrapped_coords[i]
                sites.append(pmg.Site(self.structure[i].specie, site_coord))
                site_indices.append(i)

                # Track cut bonds for attachment points
                i_frac_unwrapped = None
                for neighbor in original_undirected.neighbors(i):
                    if neighbor not in comp_indices and ((i, neighbor) in cut_edges or (neighbor, i) in cut_edges):
                        attachment_atoms.add(i)
                        cut_bonds.append((i, neighbor))

                        # Keep the peripheral functional atom that sits on the
                        # far side of the cut, so a building block is never
                        # left stripped of the chemistry at its connection
                        # point (a dioxin linker would otherwise reduce to bare
                        # benzene, a benzoxazole node to bare aldehyde
                        # carbons). The copy is merged back into a single atom
                        # when node and linker are recombined, because the
                        # assembly step de-duplicates same-element atoms that
                        # coincide in space. This mirrors the MOF convention,
                        # where carboxylate oxygens appear in both
                        # mof_nodes_lib and mof_linkers_lib.
                        #
                        # Applied to every linkage type. For imine COFs this
                        # also gives the aldehyde-derived block the partner's
                        # N (e.g. a TpAzo node becomes H3 C9 N3 O3), which
                        # shows the node's linkage environment. Assembled
                        # fragments are unaffected either way, since the
                        # duplicated atoms merge on recombination; only the
                        # standalone helper-library and OnlyNode/OnlyLinker
                        # exports differ.
                        # An azine cut is the exception: its two atoms are
                        # both nitrogen, so carrying the far one back would
                        # rebuild the N-N the cut just severed and hand the
                        # block the very -N=N-H terminus this linkage rule
                        # exists to avoid (802's node came back with three of
                        # them). One nitrogen per edge is the whole point.
                        _linkage = linkage_of.get(frozenset((i, neighbor)))
                        if _linkage is not None and _linkage[0] == 'N-N bridge':
                            continue
                        if (
                            self.structure[neighbor].specie.symbol in _CARRYOVER_ELEMENTS
                            or vinylene_partner.get(i) == neighbor
                        ):
                            if i_frac_unwrapped is None:
                                i_frac_unwrapped = self.structure.lattice.get_fractional_coords(
                                    unwrapped_coords[i]
                                )
                            _, image = self.structure.lattice.get_distance_and_image(
                                i_frac_unwrapped, self.structure[neighbor].frac_coords
                            )
                            nb_frac = self.structure[neighbor].frac_coords + image
                            carried.append(
                                (neighbor, self.structure.lattice.get_cartesian_coords(nb_frac))
                            )

            # A genuine node or linker is attached to the framework by at
            # least one severed bond. A component with none was already
            # disconnected in the parent CIF - pore solvent or a guest (929
            # carries six water/hydroxyl molecules) - and must not be exported
            # as a linker just because it has fewer than three linkages.
            if not cut_bonds:
                continue

            # Append carried heteroatoms, skipping positions already present.
            for nb_idx, nb_coord in carried:
                if any(abs(np.linalg.norm(np.asarray(s.coords) - nb_coord)) < 0.1 for s in sites):
                    continue
                sites.append(pmg.Site(self.structure[nb_idx].specie, nb_coord))
                site_indices.append(nb_idx)

            molecule = pmg.Molecule.from_sites(sites)
            
            # Create a topological graph to compute a unique hash for deduplication
            try:
                nn_strategy = JmolNN()
                mol_graph = MoleculeGraph.with_local_env_strategy(molecule, nn_strategy)
                nx_graph = mol_graph.graph.to_undirected()
                for n in nx_graph.nodes():
                    nx_graph.nodes[n]['specie'] = molecule[n].specie.symbol
                topology_hash = nx.weisfeiler_lehman_graph_hash(nx_graph, node_attr='specie')
            except Exception:
                topology_hash = ""
            
            # indices must stay aligned with the molecule's sites, which now
            # include any carried-over heteroatoms.
            sbu = SBU(molecule, indices=site_indices)
            sbu.smiles = topology_hash  # Use hash as a canonical ID for deduplication
            
            # 5. Classify as Node vs Linker by the number of distinct LINKAGES,
            # not raw attachment atoms. One chemical linkage can sever more than
            # one bond (a benzoxazole C2 loses both its C-O and its C-N bond), so
            # severed bonds that share an attachment atom or share an external
            # partner atom belong to the same linkage and must be counted once.
            # For single-bond linkages (imine, boroxine) this is identical to
            # counting attachment atoms.
            if _count_linkages(cut_bonds, linkage_of) >= 3:
                nodes.append(sbu)
            else:
                linkers.append(sbu)

        # Symmetric (k+k) COFs: when BOTH building blocks carry the same number
        # of linkages - e.g. a 3-connected triformylphenol core condensed with a
        # 3-connected tris(aminophenyl)benzene strut - the linkage count alone
        # cannot tell node from linker, so every block lands in `nodes` and
        # `linkers` comes back empty. Downstream that looks like a failed
        # partition and the caller falls through to a much cruder fallback path.
        # Break the tie by size: the node is the compact, rigid core and the
        # linker is the extended strut, so the smallest distinct block type
        # stays a node and every larger type is demoted to a linker.
        #
        # This only fires when the linkage rule produced no linkers at all AND
        # more than one distinct block type exists, so genuine single-block
        # frameworks (self-condensed boroxine) and the ordinary asymmetric case
        # (3-connected node + 2-connected linker) are untouched.
        if not linkers and len(nodes) > 1:
            groups = {}
            for sbu in nodes:
                key = (getattr(sbu, "smiles", ""), len(sbu.molecule))
                groups.setdefault(key, []).append(sbu)
            if len(groups) > 1:
                smallest = min(groups, key=lambda k: (k[1], k[0]))
                kept, demoted = [], []
                for key, members in groups.items():
                    (kept if key == smallest else demoted).extend(members)
                nodes, linkers = kept, demoted

        linkage_types = sorted({
            tag[0] if isinstance(tag, tuple) else str(tag)
            for tag in linkage_of.values()
        })
        if edges_to_remove and not linkage_types:
            # Bonds cut by the untagged rules (boroxine B-O, imine C=N).
            linkage_types = ["boroxine/imine"]
        return FragmentationResult(
            nodes, linkers,
            n_cut_bonds=len(edges_to_remove),
            linkage_types=linkage_types,
        )
