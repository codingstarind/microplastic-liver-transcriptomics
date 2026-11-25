"""Export readable question/response transcripts of the project's Claude Code sessions.

Source: the Claude Code JSONL transcripts on the analysis machine
(~/.claude/projects/C--Users-Soham-Desktop-microplastic/). For each session this writes the
user-typed prompts and Claude's written replies, in order. Tool calls, tool outputs (file contents,
command output), hidden reasoning and automatic system context are left out; the line
"[n tool calls]" marks where tools were used. The unrelated OlmoCR session (91b92e0b) is skipped,
messages shared between a forked session and its parent are written once, and email addresses are
redacted.

Run from the repository root:  python supplementary/chat_histories/export_transcripts.py
"""
import importlib.util
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("extract_prompts", HERE.parent / "prompts" / "extract_prompts.py")
ep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ep)
H = {"prompt_text": ep.prompt_text, "EMAIL": ep.EMAIL, "first_ts": ep.first_ts}
seen = set()
SRC, UNRELATED, FORKS = ep.SRC, ep.UNRELATED, ep.FORKS
for f in sorted(SRC.glob("*.jsonl"), key=lambda p: (p.stem in FORKS, H["first_ts"](p))):
    sid = f.stem
    if sid in UNRELATED:
        continue
    title, models, out, tools, first = None, Counter(), [], 0, None
    for line in open(f, encoding="utf-8"):
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        first = first or rec.get("timestamp")
        if rec.get("type") in ("custom-title", "ai-title"):
            title = rec.get("customTitle") or rec.get("aiTitle") or title
        key = rec.get("uuid")
        if key and key in seen:
            continue
        if key:
            seen.add(key)
        p = H["prompt_text"](rec)
        if p is not None:
            if tools:
                out.append(f"*[{tools} tool calls]*\n")
                tools = 0
            out.append(f"### User ({rec.get('timestamp')})\n\n{p}\n")
            continue
        if rec.get("type") != "assistant":
            continue
        msg = rec.get("message", {})
        if msg.get("model") not in (None, "<synthetic>"):
            models[msg["model"]] += 1
        for b in msg.get("content") or []:
            if not isinstance(b, dict):
                continue
            if b.get("type") == "tool_use":
                tools += 1
            elif b.get("type") == "text" and b.get("text", "").strip():
                if tools:
                    out.append(f"*[{tools} tool calls]*\n")
                    tools = 0
                out.append(f"### Claude ({msg.get('model')}, {rec.get('timestamp')})\n\n"
                           f"{H['EMAIL'].sub('[email redacted]', b['text'].strip())}\n")
    if tools:
        out.append(f"*[{tools} tool calls]*\n")
    name = HERE / f"transcript_{(first or '')[:10]}_{sid[:8]}.md"
    with open(name, "w", encoding="utf-8") as o:
        o.write(f"# Transcript: {title or 'untitled session'}\n\n- Claude Code session ID: `{sid}`\n")
        o.write(f"- Models in this export: {', '.join(f'{m} ({n} replies)' for m, n in models.most_common())}\n")
        o.write("- Prompts and written replies only; tool calls and outputs omitted.\n\n")
        o.write("\n".join(out))
    print(name.name, len(out))
