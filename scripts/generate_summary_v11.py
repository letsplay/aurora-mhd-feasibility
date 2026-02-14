"""Generate 1-page executive summary PDF — v11."""
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm, cm
from reportlab.lib.colors import HexColor, white, black
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
import os

OUTPUT = "/mnt/user-data/outputs/aurora_executive_summary_v11.pdf"
W, H = A4  # 210 x 297 mm


def draw_page(c):
    # === DARK BACKGROUND ===
    c.setFillColor(HexColor("#1a1a2e"))
    c.rect(0, 0, W, H, fill=True, stroke=False)

    # === HEADER ===
    c.setFillColor(HexColor("#00D4FF"))
    c.setFont("Helvetica-Bold", 24)
    c.drawString(2*cm, H - 2.5*cm, "AURORA-MHD")
    c.setFillColor(white)
    c.setFont("Helvetica", 11)
    c.drawString(2*cm, H - 3.2*cm,
        "Magnetohydrodynamic Thermal Protection for Starship Reentry")

    # === DIVIDER ===
    c.setStrokeColor(HexColor("#00D4FF"))
    c.setLineWidth(2)
    c.line(2*cm, H - 3.6*cm, W - 2*cm, H - 3.6*cm)

    # === PROBLEM ===
    y = H - 4.5*cm
    c.setFillColor(HexColor("#FFD700"))
    c.setFont("Helvetica-Bold", 11)
    c.drawString(2*cm, y, "PROBLEM")
    c.setFillColor(white)
    c.setFont("Helvetica", 9.5)
    y -= 0.5*cm
    c.drawString(2*cm, y,
        "Starship's ~18,000 ceramic tiles require extensive inspection after every flight.")
    y -= 0.4*cm
    c.drawString(2*cm, y,
        "A reusable magnetic heat shield could eliminate tile maintenance entirely.")

    # === APPROACH ===
    y -= 0.8*cm
    c.setFillColor(HexColor("#FFD700"))
    c.setFont("Helvetica-Bold", 11)
    c.drawString(2*cm, y, "APPROACH")
    c.setFillColor(white)
    c.setFont("Helvetica", 9.5)
    y -= 0.5*cm
    c.drawString(2*cm, y,
        "REBCO superconducting magnets generate a 2T field to deflect reentry plasma via MHD interaction.")
    y -= 0.4*cm
    c.drawString(2*cm, y,
        "Faraday power extraction from the plasma flow powers the system (self-sustaining operation).")
    y -= 0.4*cm
    c.drawString(2*cm, y,
        "Dual-model analysis: v1 (uncoupled Faraday) vs v11 (kinetic-ceiling-limited, physically corrected).")

    # === KEY RESULTS TABLE ===
    y -= 1.0*cm
    c.setFillColor(HexColor("#FFD700"))
    c.setFont("Helvetica-Bold", 11)
    c.drawString(2*cm, y, "KEY RESULTS (v11 Physics-Corrected Model)")

    y -= 0.6*cm
    col1_x = 2.2*cm
    col2_x = 9*cm
    col3_x = 13.5*cm
    row_h = 0.48*cm
    table_w = W - 4*cm

    table_data = [
        ("Metric", "v1 (naive)", "v11 (corrected)"),
        ("Coil mass (2T, 20K)", "264 kg", "~400 kg"),
        ("Flight system", "264 kg (no aux)", "~500 kg (+105 kg aux)"),
        ("Extraction model", "Uncoupled Faraday", "Kinetic ceiling-limited"),
        ("Peak extraction (LEO)", "36 MW", "~200 kW"),
        ("Energy margin", "24x (ad-hoc 1000x loss)", "~38x (physics-based)"),
        ("Self-sustaining?", "Yes (all B >= 0.5T)", "Yes (B >= 1.5T, ride-through)"),
        ("Stuart number (LEO)", "Not computed", "3,600 (S >> 1: robust deflection)"),
        ("Test suite", "58 tests", "96 tests, 100% passing"),
    ]

    # Table header
    c.setFillColor(HexColor("#00D4FF"))
    c.rect(2*cm, y - 0.1*cm, table_w, row_h, fill=True, stroke=False)
    c.setFillColor(black)
    c.setFont("Helvetica-Bold", 8.5)
    c.drawString(col1_x, y + 0.05*cm, table_data[0][0])
    c.drawString(col2_x, y + 0.05*cm, table_data[0][1])
    c.drawString(col3_x, y + 0.05*cm, table_data[0][2])

    # Table rows
    for i, (metric, v1, v11) in enumerate(table_data[1:], 1):
        row_y = y - i * row_h
        if i % 2 == 0:
            c.setFillColor(HexColor("#252545"))
            c.rect(2*cm, row_y - 0.1*cm, table_w, row_h, fill=True, stroke=False)
        c.setFillColor(white)
        c.setFont("Helvetica", 8.5)
        c.drawString(col1_x, row_y + 0.05*cm, metric)
        c.setFillColor(HexColor("#888888"))
        c.setFont("Helvetica", 8.5)
        c.drawString(col2_x, row_y + 0.05*cm, v1)
        c.setFillColor(HexColor("#00CC66"))
        c.setFont("Helvetica-Bold", 8.5)
        c.drawString(col3_x, row_y + 0.05*cm, v11)

    # === CENTRAL FINDING ===
    y -= (len(table_data)) * row_h + 0.5*cm
    c.setFillColor(HexColor("#FFD700"))
    c.setFont("Helvetica-Bold", 11)
    c.drawString(2*cm, y, "CENTRAL FINDING")
    c.setFillColor(white)
    c.setFont("Helvetica", 9.5)
    y -= 0.5*cm
    c.drawString(2*cm, y,
        "Deflection is robust (S >> 1 at all conditions). Self-sustaining extraction is feasible but marginal.")
    y -= 0.4*cm
    c.drawString(2*cm, y,
        "The v1-to-v11 progression reveals a 200x reduction in extraction power -- all physically motivated:")
    y -= 0.4*cm
    c.setFont("Helvetica", 9)
    c.drawString(2.5*cm, y,
        "Post-shock velocity (25x) + Effective area (25x) + Hall effect (1.7x) + Channel losses (3x)")
    y -= 0.4*cm
    c.setFont("Helvetica-Bold", 9.5)
    c.setFillColor(HexColor("#00CC66"))
    c.drawString(2*cm, y,
        "Coupled MHD-CFD simulation is the critical next step to resolve the remaining uncertainty.")

    # === ENVELOPE CHART ===
    y -= 0.8*cm
    chart_path = "results/figures/d3_envelope_v11.png"
    if not os.path.exists(chart_path):
        chart_path = "results/figures/d3_envelope.png"
    if os.path.exists(chart_path):
        c.setFillColor(HexColor("#FFD700"))
        c.setFont("Helvetica-Bold", 11)
        c.drawString(2*cm, y, "OPERATING ENVELOPE")
        img = ImageReader(chart_path)
        iw, ih = img.getSize()
        aspect = ih / iw
        draw_w = W - 4*cm
        draw_h = draw_w * aspect
        max_h = 5.5*cm
        if draw_h > max_h:
            draw_h = max_h
            draw_w = draw_h / aspect
        c.drawImage(chart_path,
                    2*cm, y - draw_h - 0.3*cm,
                    width=draw_w, height=draw_h)

    # === AUTHOR ===
    c.setStrokeColor(HexColor("#00D4FF"))
    c.setLineWidth(1)
    c.line(2*cm, 2.5*cm, W - 2*cm, 2.5*cm)

    c.setFillColor(HexColor("#00D4FF"))
    c.setFont("Helvetica-Bold", 10)
    c.drawString(2*cm, 1.8*cm, "Luc")
    c.setFillColor(white)
    c.setFont("Helvetica", 9)
    c.drawString(2*cm, 1.2*cm,
        "PhD Applied Physics (CEA Saclay) | 25+ years tech program management")
    c.drawString(2*cm, 0.7*cm,
        "github.com/letsplay/aurora-mhd-feasibility")

    # Version tag
    c.setFillColor(HexColor("#555555"))
    c.setFont("Helvetica", 7)
    c.drawRightString(W - 2*cm, 0.5*cm, "v2.1 — February 2026")


# === GENERATE ===
os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
c = canvas.Canvas(OUTPUT, pagesize=A4)
draw_page(c)
c.save()
print(f"Generated: {OUTPUT}")
