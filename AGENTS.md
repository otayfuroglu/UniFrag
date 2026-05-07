# Project Instructions for Codex Agents

## Required Context
Before making changes, read these root-level files:

1. `AGENTS.md`
2. `project-memory.md`
3. `project-decisions.md`
4. `project-agent-log.md`

## Memory Files
- Treat `project-memory.md` as durable project memory: architecture, workflows, commands, datasets, known risks.
- Treat `project-decisions.md` as the source of truth for durable implementation/architecture decisions.
- Treat `project-agent-log.md` as the chronological handoff log between agents. It is historical, not policy.

## Workflow
- Read `project-memory.md`, `project-decisions.md`, and recent `project-agent-log.md` entries before editing.
- Follow existing architecture and conventions.
- Do not introduce new frameworks without a clear reason.
- Keep code consistent with existing style.
- After finishing work, update `project-agent-log.md` with what changed, files touched, validation, and follow-up risks.
- When durable project context changes, update `project-memory.md`.
- When a durable implementation decision is made or changed, update `project-decisions.md`.

## Conflict Rule
If instructions conflict:
1. `AGENTS.md` wins.
2. `project-decisions.md` wins next.
3. `project-memory.md` wins next.
4. `project-agent-log.md` is historical only.
