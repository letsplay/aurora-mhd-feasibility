"""Generate 1-page executive summary PDF."""
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm, cm
from reportlab.lib.colors import HexColor, white, black
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
import os

OUTPUT = "results/figures/aurora_executive_summary.pdf"
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
    c.setFont("Helvetica", 12)
    c.drawString(2*cm, H - 3.2*cm,
        "Can MHD Thermal Protection Be Self-Sustaining for Starship Reentry?")

    # === DIVIDER ===
    c.setStrokeColor(HexColor("#00D4FF"))
    c.setLineWidth(2)
    c.line(2*cm, H - 3.6*cm, W - 2*cm, H - 3.6*cm)

    # === PROBLEM (2 sentences) ===
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

    # === APPROACH (3 sentences) ===
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
        "Faraday power extraction from the plasma flow powers the system (self-sustaining).")
    y -= 0.4*cm
    c.drawString(2*cm, y,
        "Physics-based models validated against Biot-Savart, NRLMSISE-00, Rankine-Hugoniot, and Saha equation.")

    # === KEY RESULTS TABLE ===
    y -= 1.0*cm
    c.setFillColor(HexColor("#FFD700"))
    c.setFont("Helvetica-Bold", 11)
    c.drawString(2*cm, y, "KEY RESULTS")

    y -= 0.6*cm
    table_data = [
        ("Metric", "Value"),
        ("Magnet mass (2T, 20K)", "264 kg"),
        ("Energy closure", "Self-sustaining at all B >= 0.5T"),
        ("Safety margin (1000x loss)", "24x (LEO) to 77x (Mars)"),
        ("MHD effective zone", "v > 3 km/s, h > 40 km"),
        ("Test suite", "58 tests, 100% passing"),
    ]

    # Table header
    c.setFillColor(HexColor("#00D4FF"))
    c.rect(2*cm, y - 0.1*cm, W - 4*cm, 0.5*cm, fill=True, stroke=False)
    c.setFillColor(black)
    c.setFont("Helvetica-Bold", 9)
    c.drawString(2.2*cm, y + 0.05*cm, table_data[0][0])
    c.drawString(9*cm, y + 0.05*cm, table_data[0][1])

    # Table rows
    for i, (metric, value) in enumerate(table_data[1:], 1):
        row_y = y - i * 0.5*cm
        if i % 2 == 0:
            c.setFillColor(HexColor("#252545"))
            c.rect(2*cm, row_y - 0.1*cm, W - 4*cm, 0.5*cm, fill=True, stroke=False)
        c.setFillColor(white)
        c.setFont("Helvetica", 9)
        c.drawString(2.2*cm, row_y + 0.05*cm, metric)
        c.setFont("Helvetica-Bold", 9)
        c.drawString(9*cm, row_y + 0.05*cm, value)

    # === ENVELOPE CHART ===
    chart_y = y - 4.0*cm
    chart_path = "results/figures/d3_envelope.png"
    if os.path.exists(chart_path):
        c.setFillColor(HexColor("#FFD700"))
        c.setFont("Helvetica-Bold", 11)
        c.drawString(2*cm, chart_y, "OPERATING ENVELOPE")
        img = ImageReader(chart_path)
        iw, ih = img.getSize()
        aspect = ih / iw
        draw_w = W - 4*cm
        draw_h = draw_w * aspect
        if draw_h > 7*cm:
            draw_h = 7*cm
            draw_w = draw_h / aspect
        c.drawImage(chart_path,
                    2*cm, chart_y - draw_h - 0.3*cm,
                    width=draw_w, height=draw_h)
        bottom_y = chart_y - draw_h - 0.5*cm
    else:
        bottom_y = chart_y - 1*cm
        c.setFillColor(HexColor("#FF3344"))
        c.setFont("Helvetica", 9)
        c.drawString(2*cm, bottom_y, "[Chart not found — run notebooks first]")

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

# === GENERATE ===
os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
c = canvas.Canvas(OUTPUT, pagesize=A4)
draw_page(c)
c.save()
print(f"Generated: {OUTPUT}")