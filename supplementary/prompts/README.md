# Prompts

The prompts the team typed to Claude Code during the project are exported here, one file per session, in time order with UTC time stamps. Pasted instructions are included as they were pasted.

The exports leave out:
- automatic notices (sub-agent completion messages, CI notices);
- Claude's replies, which are in [`../chat_histories/`](../chat_histories/).

Email addresses are redacted.

| File | Prompts | Topic |
|---|---|---|
| [session_2026-11-02_949232bd.md](session_2026-11-02_949232bd.md) | 24 | Set-up, Phase 1, methodology, re-quantification, report draft |
| [session_2026-11-02_03c789f5.md](session_2026-11-02_03c789f5.md) | 3 | Phase 2 factorial modelling and pathways |
| [session_2026-11-20_9228d9bf.md](session_2026-11-20_9228d9bf.md) | 8 | Literature review |
| [session_2026-11-21_eb0c5c48.md](session_2026-11-21_eb0c5c48.md) | 1 | Published-gene replication |
| [session_2026-11-21_5dbfd995.md](session_2026-11-21_5dbfd995.md) | 3 | Research strategy |
| [session_2026-11-22_59dc8d15.md](session_2026-11-22_59dc8d15.md) | 39 | Phases 2–5C, follow-ups, figures, repository organisation |
| [session_2026-11-22_094de3c2.md](session_2026-11-22_094de3c2.md) | 1 | Fork for the Section 3.1 figures |

`prompt_index.json` lists the same sessions as machine-readable data. `extract_prompts.py` rebuilds this folder from the local transcripts.
