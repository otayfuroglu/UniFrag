import os
from pathlib import Path
from pathlib import Path
import pymatgen.core as pmg
from pymatgen.analysis.graphs import StructureGraph, MoleculeGraph
from pymatgen.analysis.local_env import JmolNN
import networkx as nx

class SBU:
    """Represents a Secondary Building Unit (Node or Linker) in a COF."""
    def __init__(self, molecule, indices=None):
        self.molecule = molecule # pymatgen.core.Molecule
        self.indices = indices
        self.smiles = "" # Placeholder for SMILES, in practice calculated via RDKit/OpenBabel
        self.molar_mass = molecule.composition.weight if molecule else 0

class FragmentationResult:
    """Holds the resulting nodes and linkers from COF fragmentation."""
    def __init__(self, nodes, linkers):
        self.nodes = nodes
        self.linkers = linkers

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
        # Create an undirected version for component analysis
        undirected_graph = sg.graph.to_undirected()
        
        # Create a copy to preserve original connections for neighbor analysis
        original_undirected = undirected_graph.copy()
        
        # Helper: count heavy (non-H) neighbors of a node
        def heavy_degree(node):
            return sum(1 for nb in undirected_graph.neighbors(node)
                       if self.structure[nb].specie.symbol != 'H')

        # Helper: check if an edge is in a small ring (e.g. porphyrin pyrrole ring)
        def in_small_ring(u, v, max_len=6):
            G_temp = undirected_graph.copy()
            G_temp.remove_edge(u, v)
            try:
                path = nx.shortest_path(G_temp, u, v)
                return len(path) <= max_len
            except nx.NetworkXNoPath:
                return False

        for u, v, data in undirected_graph.edges(data=True):
            atom_u = self.structure[u].specie.symbol
            atom_v = self.structure[v].specie.symbol
            
            bond_pair = set([atom_u, atom_v])
            
            # Heuristic: break common COF forming bonds.
            if bond_pair == {'B', 'O'}:
                edges_to_remove.append((u, v))
            elif bond_pair == {'C', 'N'}:
                # Cut on the aromatic C side (degree >= 3) so that N stays with the linker.
                # E.g. in imine COFs: BDA-CH=N-Ar(node) → cut Ar(node)-N bond.
                # Protect C-N bonds that are part of small rings (porphyrin, triazine, pyridine).
                c_idx = u if atom_u == 'C' else v
                if heavy_degree(c_idx) >= 3 and not in_small_ring(u, v):
                    edges_to_remove.append((u, v))
                
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
            
            for i in comp_indices:
                # Use the unwrapped coordinates instead of raw unit cell coordinates
                site_coord = unwrapped_coords[i]
                sites.append(pmg.Site(self.structure[i].specie, site_coord))
                
                # Track cut bonds for attachment points
                for neighbor in original_undirected.neighbors(i):
                    if neighbor not in comp_indices and ((i, neighbor) in cut_edges or (neighbor, i) in cut_edges):
                        attachment_atoms.add(i)
                            
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
            
            sbu = SBU(molecule, indices=comp_indices)
            sbu.smiles = topology_hash  # Use hash as a canonical ID for deduplication
            
            # 5. Classify as Node vs Linker based on number of attachment atoms
            if len(attachment_atoms) >= 3:
                nodes.append(sbu)
            else:
                linkers.append(sbu)
                
        return FragmentationResult(nodes, linkers)
