# ChemSmart Agent Architecture Audit

This is the measured pre-refactor baseline. It is not evidence that the structural completion gates have been met.

- Branch: `codex/agent-codebase-simplification`
- HEAD: `95faa53cf1b231f96f0f075a977f46be1c6e1ebe`

| Metric | Value |
|---|---:|
| Python files | 151 |
| Physical lines | 44691 |
| Functions and methods | 1694 |
| Classes | 202 |
| Files over 800 lines | 14 |
| Classes over 500 lines | 4 |
| Functions over 100 lines | 39 |
| Import cycles | 1 |
| Ruff default violations | 0 |
| Ruff C901 violations | 65 |
| Pylint duplicate groups | 8 |
| Branch-aware test coverage | 74.63% |

## Oversized Files

| File | Lines |
|---|---:|
| `chemsmart/agent/tui/screens/chat.py` | 4628 |
| `chemsmart/agent/core.py` | 2621 |
| `chemsmart/agent/tools.py` | 1633 |
| `chemsmart/agent/synthesis.py` | 1551 |
| `chemsmart/agent/project_yaml.py` | 1256 |
| `chemsmart/agent/runtime/calculations.py` | 1136 |
| `chemsmart/agent/cli.py` | 962 |
| `chemsmart/agent/loop.py` | 950 |
| `chemsmart/agent/services/conversation_memory.py` | 947 |
| `chemsmart/agent/model_command_parser.py` | 927 |
| `chemsmart/agent/harness/command_semantics.py` | 904 |
| `chemsmart/agent/harness/generated_invariants.py` | 883 |
| `chemsmart/agent/tools_command.py` | 831 |
| `chemsmart/agent/harness/command_contracts.py` | 818 |

## Long Classes

| Class | File | Lines |
|---|---|---:|
| `ChatScreen` | `chemsmart/agent/tui/screens/chat.py` | 3961 |
| `AgentSession` | `chemsmart/agent/core.py` | 1776 |
| `SynthesisSession` | `chemsmart/agent/synthesis.py` | 862 |
| `ToolLoop` | `chemsmart/agent/loop.py` | 752 |

## Longest Functions

| Function | File | Lines |
|---|---|---:|
| `AgentSession.run_loop` | `chemsmart/agent/core.py` | 357 |
| `ToolLoop.run_turn` | `chemsmart/agent/loop.py` | 313 |
| `evaluate_command_semantics` | `chemsmart/agent/harness/command_semantics.py` | 240 |
| `_summarize_tool_use_result` | `chemsmart/agent/services/conversation_memory.py` | 235 |
| `AgentSession._continue_run` | `chemsmart/agent/core.py` | 227 |
| `ChatScreen._apply_log_entry` | `chemsmart/agent/tui/screens/chat.py` | 215 |
| `parse_decision_event` | `chemsmart/agent/tui/events.py` | 207 |
| `ToolRegistry.default` | `chemsmart/agent/registry.py` | 188 |
| `parse_model_command` | `chemsmart/agent/model_command_parser.py` | 183 |
| `ChatScreen._handle_slash_command` | `chemsmart/agent/tui/screens/chat.py` | 182 |
| `execute_observed_process` | `chemsmart/agent/runtime/calculations.py` | 179 |
| `_qmmm_contract_issues` | `chemsmart/agent/harness/command_contracts.py` | 175 |
| `ChatScreen._publish_synthesis_result` | `chemsmart/agent/tui/screens/chat.py` | 173 |
| `check_command_contracts` | `chemsmart/agent/harness/command_contracts.py` | 170 |
| `_execution_terminal_state` | `chemsmart/agent/tools_command.py` | 166 |
| `AgentSession._finalize_session` | `chemsmart/agent/core.py` | 146 |
| `inspect_output` | `chemsmart/agent/runtime/calculations.py` | 146 |
| `ChatScreen._apply_tool_use_event` | `chemsmart/agent/tui/screens/chat.py` | 144 |
| `ChatScreen.run_slash_tool_request` | `chemsmart/agent/tui/screens/chat.py` | 139 |
| `build_sub_intent_assertions` | `chemsmart/agent/harness/sub_intent.py` | 135 |

## Dependency Hotspots

| Module | Fan-in | Fan-out |
|---|---:|---:|
| `chemsmart.agent.tui.screens.chat` | 1 | 33 |
| `chemsmart.agent.core` | 8 | 17 |
| `chemsmart.agent.tui.widgets.cells` | 1 | 19 |
| `chemsmart.agent.tui.widgets.cells.base` | 19 | 0 |
| `chemsmart.agent.synthesis` | 3 | 15 |
| `chemsmart.agent.wizard` | 1 | 16 |
| `chemsmart.agent.cli` | 1 | 12 |
| `chemsmart.agent.tools_command` | 1 | 11 |
| `chemsmart.agent.wizard.orchestrator` | 2 | 9 |
| `chemsmart.agent.wizard.probe` | 9 | 2 |
| `chemsmart.agent.providers` | 4 | 6 |
| `chemsmart.agent.wizard.refresh` | 3 | 7 |
| `chemsmart.agent.harness.workflow_state` | 9 | 0 |
| `chemsmart.agent.runtime.orchestrator` | 2 | 7 |
| `chemsmart.agent.tui.widgets.popups` | 1 | 8 |

## Import Cycles

- `chemsmart.agent.cli -> chemsmart.agent.tui -> chemsmart.agent.tui.app -> chemsmart.agent.tui.screens.chat`

## Duplicate Groups

1. `chemsmart.agent.cli:[667:687]` and `chemsmart.agent.tui.screens.chat:[3387:3407]`
2. `chemsmart.agent.local.adapter:[181:199]` and `chemsmart.agent.synthesis:[1150:1168]`
3. `chemsmart.agent.wizard.project:[155:177]` and `chemsmart.agent.wizard.scratch:[145:167]`
4. `chemsmart.agent.harness.intent:[393:412]` and `chemsmart.agent.harness.sub_intent:[74:107]`
5. `chemsmart.agent.wizard.project:[139:155]` and `chemsmart.agent.wizard.software:[407:427]`
6. `chemsmart.agent.harness.intent:[42:53]` and `chemsmart.agent.tools_command:[331:342]`
7. `chemsmart.agent.wizard.scratch:[153:167]` and `chemsmart.agent.wizard.software:[548:568]`
8. `chemsmart.agent.harness.intent:[314:324]` and `chemsmart.agent.tools_command:[331:341]`

## Completion Targets

The final branch must reduce every production agent file to 800 lines or fewer, every class to 500 lines or fewer, and every function to 100 lines or fewer. C901 findings, duplicate groups, and import cycles must all reach zero without increasing total agent production LOC or changing the frozen behavior contracts.
