"""Gradio UI skeleton — upload an Embark CSV, get the two-layer report.

Runs against WILDTYPE_MODE from .env (mock by default, so this is fully
clickable tonight with no sponsor keys). Deploy target tomorrow: Modal, per
the proposal's hosting plan — `modal deploy` wraps this same app, no UI
rewrite needed, just an entrypoint file (add tomorrow once Modal access
is live).

    python3 -m wildtype.ui.app
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from dotenv import load_dotenv

load_dotenv()

import gradio as gr

from wildtype.agent.loop import maybe_escalate_to_af3, render_report, run_pipeline
from wildtype.parsers.embark_csv import parse_report
from wildtype.tools.base import get_toolset

_RISK_LABEL = {"at_risk": "⚠️ At risk", "carrier": "🧬 Carrier", "clear": "✅ Clear", "unknown": "? Unknown"}


def analyze(csv_file, breed_name: str, species: str):
    mode = os.environ.get("WILDTYPE_MODE", "mock")
    toolset = get_toolset(mode)

    if csv_file is not None:
        report = parse_report(csv_file.name)
    elif breed_name.strip():
        # Breed-name-only mode is a nice-to-have (proposal: needs Paperclip
        # to fetch known-disease variants for the breed). Until that's
        # wired, surface the gap honestly instead of pretending to run it.
        return (
            f"Breed-name-only mode (\"{breed_name}\") isn't wired to Paperclip yet — "
            "nice-to-have per the plan. Upload an Embark CSV for now.",
            "",
        )
    else:
        return "Upload an Embark CSV or enter a breed name.", ""

    analyses = run_pipeline(report, toolset, also_include_genes=["ABCB1"])
    if not analyses:
        return f"No flagged variants found for {report.name}.", ""

    plain_sections, science_sections = [], []
    for analysis in analyses:
        analysis = maybe_escalate_to_af3(analysis, toolset)
        full_report = render_report(analysis)
        label = _RISK_LABEL.get(analysis.variant.risk_level, analysis.variant.risk_level)
        header = f"### {label} — {analysis.variant.gene_hint}\n{analysis.variant.name}\n"
        if "Science layer" in full_report:
            plain, _, science = full_report.partition("Science layer")
            plain_sections.append(header + plain.replace("Plain layer", "").strip())
            science_sections.append(header + "\n" + science.strip())
        else:
            plain_sections.append(header + full_report)

    subject = f"## {report.name} — {report.breed_name} ({species})\n\n"
    return subject + "\n\n---\n\n".join(plain_sections), "\n\n---\n\n".join(science_sections)


with gr.Blocks(title="WildType") as demo:
    gr.Markdown("# WildType\nSpecies-agnostic genomic variant interpretation.")

    with gr.Row():
        csv_input = gr.File(label="Upload Embark CSV", file_types=[".csv"])
        breed_input = gr.Textbox(label="…or enter breed name (nice-to-have, not wired yet)")
    species_input = gr.Radio(["Dog", "Cat", "Other"], value="Dog", label="Species")
    run_btn = gr.Button("Run Analysis", variant="primary")

    gr.Markdown("## Results")
    with gr.Row():
        plain_output = gr.Markdown(label="Pet Owner View")
        science_output = gr.Markdown(label="Researcher View")

    run_btn.click(analyze, inputs=[csv_input, breed_input, species_input], outputs=[plain_output, science_output])

if __name__ == "__main__":
    demo.launch()
