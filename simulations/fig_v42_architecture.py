"""
V42 Model Architecture Diagram — MYCN-EZH2-CyclinD1-HH Cell Cycle Model.
Generates a publication-quality schematic of all model components and connections.
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

fig, ax = plt.subplots(figsize=(16, 11))
ax.set_xlim(0, 16)
ax.set_ylim(0, 11)
ax.set_aspect('equal')
ax.axis('off')

# =====================================================================
# COLOR SCHEME
# =====================================================================
C_HH = '#d4edda'        # light green (HH module bg)
C_HH_BORDER = '#28a745'
C_MYCN = '#fff3cd'       # light yellow/orange (MYCN module bg)
C_MYCN_BORDER = '#e67e22'
C_EZH2 = '#e8daef'       # light purple (EZH2 module bg)
C_EZH2_BORDER = '#8e44ad'
C_CC = '#d6eaf8'         # light blue (cell cycle bg)
C_CC_BORDER = '#2980b9'
C_DRUG = '#e74c3c'       # red for drugs
C_ACTIVATE = '#27ae60'   # green arrows (activation)
C_INHIBIT = '#e74c3c'    # red arrows (inhibition)
C_TXLATE = '#2c3e50'     # dark (transcription/translation)
C_FEEDBACK = '#8e44ad'   # purple (feedback)
C_NODE = '#ffffff'        # white node fill

# =====================================================================
# HELPER FUNCTIONS
# =====================================================================

def draw_node(x, y, text, color='#ecf0f1', border='#2c3e50', fontsize=9,
              width=1.6, height=0.55, bold=False, fontstyle='normal'):
    """Draw a rounded rectangle node with text."""
    box = FancyBboxPatch((x - width/2, y - height/2), width, height,
                         boxstyle="round,pad=0.08", facecolor=color,
                         edgecolor=border, linewidth=1.5, zorder=3)
    ax.add_patch(box)
    weight = 'bold' if bold else 'normal'
    ax.text(x, y, text, ha='center', va='center', fontsize=fontsize,
            fontweight=weight, fontstyle=fontstyle, zorder=4)
    return (x, y, width, height)


def draw_module_bg(x, y, w, h, color, border, label, label_pos='top'):
    """Draw a module background rectangle with label."""
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15",
                         facecolor=color, edgecolor=border,
                         linewidth=2, alpha=0.4, zorder=1)
    ax.add_patch(box)
    if label_pos == 'top':
        ax.text(x + w/2, y + h + 0.15, label, ha='center', va='bottom',
                fontsize=11, fontweight='bold', color=border, zorder=2)
    elif label_pos == 'bottom':
        ax.text(x + w/2, y - 0.15, label, ha='center', va='top',
                fontsize=11, fontweight='bold', color=border, zorder=2)


def arrow(x1, y1, x2, y2, color=C_ACTIVATE, style='->', lw=1.5,
          connectionstyle='arc3,rad=0', linestyle='-', shrinkA=12, shrinkB=12):
    """Draw an arrow between two points."""
    if style == '-|':
        # Inhibition: flat head
        arr = FancyArrowPatch((x1, y1), (x2, y2),
                              arrowstyle='-[', mutation_scale=8,
                              connectionstyle=connectionstyle,
                              color=color, linewidth=lw, linestyle=linestyle,
                              shrinkA=shrinkA, shrinkB=shrinkB, zorder=2)
    else:
        arr = FancyArrowPatch((x1, y1), (x2, y2),
                              arrowstyle='->', mutation_scale=15,
                              connectionstyle=connectionstyle,
                              color=color, linewidth=lw, linestyle=linestyle,
                              shrinkA=shrinkA, shrinkB=shrinkB, zorder=2)
    ax.add_patch(arr)


def inhibit_arrow(x1, y1, x2, y2, color=C_INHIBIT, lw=1.8,
                  connectionstyle='arc3,rad=0', linestyle='-'):
    """Draw an inhibition arrow (flat bar head)."""
    arrow(x1, y1, x2, y2, color=color, style='-|', lw=lw,
          connectionstyle=connectionstyle, linestyle=linestyle)


# =====================================================================
# MODULE BACKGROUNDS
# =====================================================================

# Hedgehog Signaling module
draw_module_bg(0.3, 6.0, 5.8, 4.2, C_HH, C_HH_BORDER,
               'Hedgehog Signaling')

# MYCN module
draw_module_bg(6.6, 6.8, 3.2, 3.4, C_MYCN, C_MYCN_BORDER,
               'MYCN Module')

# EZH2 module
draw_module_bg(10.3, 6.8, 3.2, 3.4, C_EZH2, C_EZH2_BORDER,
               'EZH2 Module')

# CycD1 integration zone
draw_module_bg(3.5, 5.1, 7.0, 1.2, '#fdebd0', '#e67e22',
               'CycD1 Integration', label_pos='bottom')

# Cell Cycle Engine
draw_module_bg(0.3, 0.5, 15.2, 4.2, C_CC, C_CC_BORDER,
               'Cell Cycle Engine (Gerard 2009)', label_pos='bottom')

# =====================================================================
# NODES — HEDGEHOG PATHWAY
# =====================================================================

# SHH ligand
draw_node(1.5, 9.6, 'SHH', color='#a9dfbf', border=C_HH_BORDER, fontsize=10, bold=True)

# Ptch1
draw_node(1.5, 8.5, 'Ptch1', color=C_NODE, border=C_HH_BORDER, fontsize=9)

# Smo
draw_node(3.5, 8.5, 'Smo', color=C_NODE, border=C_HH_BORDER, fontsize=9, bold=True)

# Gli processing
draw_node(3.5, 7.2, 'Gli_act', color='#a9dfbf', border=C_HH_BORDER, fontsize=8.5)
draw_node(1.5, 7.2, 'Gli_rep', color='#f5b7b1', border='#c0392b', fontsize=8.5)

# Gli1
draw_node(5.0, 7.2, 'Gli1', color='#a9dfbf', border=C_HH_BORDER, fontsize=9)

# =====================================================================
# NODES — MYCN MODULE
# =====================================================================

draw_node(8.2, 9.0, 'MYCN', color='#fdebd0', border=C_MYCN_BORDER,
          fontsize=10, bold=True, width=1.8)
ax.text(8.2, 8.35, 'amplified in MB\n(Hill n=3)', ha='center', va='center',
        fontsize=7.5, fontstyle='italic', color='#7f6c5c', zorder=4)

# =====================================================================
# NODES — EZH2 MODULE
# =====================================================================

draw_node(11.9, 9.0, 'EZH2', color='#d7bde2', border=C_EZH2_BORDER,
          fontsize=10, bold=True, width=1.8)
ax.text(11.9, 8.35, 'H3K27me3\nrepresses CycD1', ha='center', va='center',
        fontsize=7.5, fontstyle='italic', color='#6c3483', zorder=4)

# =====================================================================
# NODES — DRUG TARGETS
# =====================================================================

draw_node(5.0, 9.0, 'GDC0449', color='#fadbd8', border=C_DRUG,
          fontsize=8.5, bold=True, width=1.7)
ax.text(5.0, 9.45, 'SMO inhibitor', ha='center', va='bottom',
        fontsize=7, color=C_DRUG, fontstyle='italic')

draw_node(14.3, 9.0, 'EZH2i', color='#fadbd8', border=C_DRUG,
          fontsize=8.5, bold=True, width=1.5)
ax.text(14.3, 9.45, 'PRC2 inhibitor', ha='center', va='bottom',
        fontsize=7, color=C_DRUG, fontstyle='italic')

# =====================================================================
# NODE — CycD1 mRNA (INTEGRATION NODE)
# =====================================================================

draw_node(7.0, 5.7, 'CycD1 mRNA', color='#fdebd0', border='#d35400',
          fontsize=10, bold=True, width=2.2, height=0.65)
ax.text(7.0, 5.15, '= (basal + Gli + MYCN)\n    x EZH2 repression',
        ha='center', va='center', fontsize=7.5, fontstyle='italic',
        color='#784212', zorder=4,
        bbox=dict(boxstyle='round,pad=0.2', facecolor='#fef9e7', alpha=0.7))

# =====================================================================
# NODES — CELL CYCLE ENGINE
# =====================================================================

# G1 phase
draw_node(2.0, 3.6, 'CycD1\nprotein', color=C_NODE, border=C_CC_BORDER, fontsize=8, height=0.6)
draw_node(4.0, 3.6, 'CycD/\nCDK4', color='#aed6f1', border=C_CC_BORDER, fontsize=8.5, bold=True, height=0.6)

# Rb-E2F switch
draw_node(6.5, 3.6, 'pRB/E2F\nswitch', color='#aed6f1', border=C_CC_BORDER,
          fontsize=8.5, bold=True, height=0.6, width=1.8)

# S phase
draw_node(9.0, 3.6, 'CycE/\nCDK2', color='#aed6f1', border=C_CC_BORDER, fontsize=8.5, bold=True, height=0.6)

# G2 phase
draw_node(11.0, 3.6, 'CycA/\nCDK2', color='#aed6f1', border=C_CC_BORDER, fontsize=8.5, bold=True, height=0.6)

# M phase
draw_node(13.0, 3.6, 'CycB/\nCDK1', color='#aed6f1', border=C_CC_BORDER, fontsize=8.5, bold=True, height=0.6)

# Division output
draw_node(14.8, 3.6, 'Division', color='#abebc6', border='#1e8449',
          fontsize=9, bold=True, width=1.5)

# p27 barrier
draw_node(7.5, 1.5, 'p27/Kip1', color='#f5b7b1', border='#c0392b',
          fontsize=9, bold=True, width=1.8)
ax.text(7.5, 0.9, 'CDK inhibitor\n(G0 arrest barrier)', ha='center', va='center',
        fontsize=7.5, fontstyle='italic', color='#922b21', zorder=4)

# APC/Cdh1 and Skp2
draw_node(3.5, 1.5, 'Cdh1\n(APC/C)', color=C_NODE, border=C_CC_BORDER, fontsize=7.5, height=0.55)
draw_node(5.5, 1.5, 'Skp2', color=C_NODE, border=C_CC_BORDER, fontsize=8)

# Cdc20
draw_node(11.5, 1.5, 'Cdc20\n(APC/C)', color=C_NODE, border=C_CC_BORDER, fontsize=7.5, height=0.55)

# Cdc25 phosphatases
draw_node(9.0, 2.3, 'Cdc25', color=C_NODE, border=C_CC_BORDER, fontsize=7.5, width=1.2, height=0.4)

# E2F (free)
draw_node(6.5, 2.5, 'E2F\n(free)', color='#d5f5e3', border=C_CC_BORDER, fontsize=8, height=0.55)

# Phase labels
ax.text(2.0, 4.35, 'G1', ha='center', fontsize=8, color='#5d6d7e', fontstyle='italic')
ax.text(4.0, 4.35, 'G1', ha='center', fontsize=8, color='#5d6d7e', fontstyle='italic')
ax.text(9.0, 4.35, 'G1/S', ha='center', fontsize=8, color='#5d6d7e', fontstyle='italic')
ax.text(11.0, 4.35, 'S/G2', ha='center', fontsize=8, color='#5d6d7e', fontstyle='italic')
ax.text(13.0, 4.35, 'M', ha='center', fontsize=8, color='#5d6d7e', fontstyle='italic')

# =====================================================================
# ARROWS — HEDGEHOG PATHWAY
# =====================================================================

# SHH → sequesters Ptch1
arrow(1.5, 9.3, 1.5, 8.8, color=C_ACTIVATE)

# Ptch1 ⊣ Smo
inhibit_arrow(2.3, 8.5, 2.7, 8.5, color=C_INHIBIT)

# Smo → Gli_act (activation)
arrow(3.5, 8.2, 3.5, 7.5, color=C_ACTIVATE)

# Gli_act ↔ Gli_rep (interconversion)
arrow(3.0, 7.2, 2.1, 7.2, color='#7f8c8d', lw=1.0)
arrow(2.0, 7.05, 2.9, 7.05, color='#7f8c8d', lw=1.0)
ax.text(2.5, 7.45, 'Smo-dependent\nswitch', ha='center', fontsize=6.5,
        color='#5d6d7e', fontstyle='italic')

# Gli_act → Gli1 (positive feedback)
arrow(4.0, 7.2, 4.4, 7.2, color=C_ACTIVATE, lw=1.3)

# GDC0449 ⊣ Smo
inhibit_arrow(4.6, 8.8, 4.0, 8.6, color=C_DRUG, lw=2.0)

# =====================================================================
# ARROWS — MYCN
# =====================================================================

# Gli → MYCN (small component)
arrow(5.5, 7.4, 7.3, 8.8, color=C_ACTIVATE, lw=1.2, linestyle='--',
      connectionstyle='arc3,rad=-0.15')
ax.text(6.0, 8.3, 'Gli-dependent\n(~34%)', fontsize=6.5, color='#7f6c5c',
        fontstyle='italic', ha='center', rotation=30)

# MYCN autonomous (label)
ax.text(8.2, 9.75, 'autonomous\n(GDC-resistant)', ha='center', va='bottom',
        fontsize=7, fontstyle='italic', color=C_MYCN_BORDER)

# MYCN → CycD1 mRNA
arrow(8.2, 8.0, 7.6, 6.05, color=C_MYCN_BORDER, lw=2.0,
      connectionstyle='arc3,rad=0.1')

# =====================================================================
# ARROWS — TO CycD1 mRNA
# =====================================================================

# Gli_act/Gli1 → CycD1 mRNA
arrow(4.2, 6.9, 5.9, 5.9, color=C_ACTIVATE, lw=1.8,
      connectionstyle='arc3,rad=-0.15')

# EZH2 ⊣ CycD1 mRNA
inhibit_arrow(11.0, 8.8, 8.2, 6.0, color=C_EZH2_BORDER, lw=2.0,
              connectionstyle='arc3,rad=-0.15')

# EZH2i ⊣ EZH2
inhibit_arrow(13.5, 9.0, 12.8, 9.0, color=C_DRUG, lw=2.0)

# =====================================================================
# ARROWS — CycD1 mRNA → CELL CYCLE
# =====================================================================

# CycD1 mRNA → CycD1 protein (translation)
arrow(5.9, 5.5, 2.5, 3.95, color=C_TXLATE, lw=1.5,
      connectionstyle='arc3,rad=0.1')
ax.text(3.8, 4.85, 'translation', fontsize=7, color=C_TXLATE, fontstyle='italic',
        rotation=18)

# CycD1 → CycD/CDK4
arrow(2.8, 3.6, 3.2, 3.6, color=C_ACTIVATE, lw=1.5)

# CycD/CDK4 → pRB/E2F
arrow(4.8, 3.6, 5.5, 3.6, color=C_ACTIVATE, lw=1.8)
ax.text(5.1, 3.95, 'pRB\nphosph.', fontsize=6.5, ha='center',
        color=C_CC_BORDER, fontstyle='italic')

# pRB/E2F → CycE/CDK2
arrow(7.4, 3.6, 8.2, 3.6, color=C_ACTIVATE, lw=1.8)
ax.text(7.8, 3.95, 'E2F\nrelease', fontsize=6.5, ha='center',
        color=C_CC_BORDER, fontstyle='italic')

# CycE/CDK2 → CycA/CDK2
arrow(9.8, 3.6, 10.2, 3.6, color=C_ACTIVATE, lw=1.5)

# CycA/CDK2 → CycB/CDK1
arrow(11.8, 3.6, 12.2, 3.6, color=C_ACTIVATE, lw=1.5)

# CycB/CDK1 → Division
arrow(13.8, 3.6, 14.0, 3.6, color=C_ACTIVATE, lw=1.8)

# =====================================================================
# ARROWS — p27 INHIBITION
# =====================================================================

# p27 ⊣ CycD/CDK4
inhibit_arrow(6.6, 1.7, 4.5, 3.3, color=C_INHIBIT, lw=1.5,
              connectionstyle='arc3,rad=0.1')

# p27 ⊣ CycE/CDK2
inhibit_arrow(7.8, 1.8, 8.8, 3.3, color=C_INHIBIT, lw=1.5,
              connectionstyle='arc3,rad=-0.05')

# p27 ⊣ CycA/CDK2
inhibit_arrow(8.3, 1.7, 10.7, 3.3, color=C_INHIBIT, lw=1.2, linestyle='--',
              connectionstyle='arc3,rad=-0.1')

# =====================================================================
# ARROWS — FEEDBACK LOOPS
# =====================================================================

# E2F → EZH2 (positive regulation)
arrow(7.0, 2.8, 11.2, 7.5, color=C_FEEDBACK, lw=1.5, linestyle='--',
      connectionstyle='arc3,rad=-0.3')
ax.text(10.5, 5.0, 'E2F / pRBpp\ndrive EZH2', fontsize=7, color=C_FEEDBACK,
        fontstyle='italic', ha='center', rotation=55)

# E2F → CycE synthesis
arrow(6.8, 2.8, 8.5, 3.3, color=C_CC_BORDER, lw=1.2, linestyle='--',
      connectionstyle='arc3,rad=-0.05')

# Skp2 → degrades p27
arrow(5.8, 1.3, 6.6, 1.3, color=C_ACTIVATE, lw=1.2)
ax.text(6.2, 1.05, 'degrades', fontsize=6.5, ha='center', color='#5d6d7e')

# Cdh1 ⊣ Skp2
inhibit_arrow(4.2, 1.5, 4.8, 1.5, color=C_INHIBIT, lw=1.2)

# CycB/CDK1 → Cdc20 activation
arrow(13.0, 3.3, 12.0, 1.8, color=C_ACTIVATE, lw=1.2, linestyle='--',
      connectionstyle='arc3,rad=-0.1')

# Cdc20 → degrades CycB (negative feedback)
arrow(11.2, 1.8, 12.5, 3.3, color=C_INHIBIT, lw=1.2, linestyle='--',
      connectionstyle='arc3,rad=-0.2')

# Cdc25 activates CDKs
arrow(9.0, 2.5, 9.0, 3.3, color=C_ACTIVATE, lw=1.0, linestyle='--')

# =====================================================================
# TITLE AND LEGEND
# =====================================================================

ax.text(8.0, 10.85, 'V42 Model Architecture — MYCN-EZH2-CyclinD1-Hedgehog Cell Cycle',
        ha='center', va='center', fontsize=14, fontweight='bold', color='#2c3e50')

# Legend
leg_x, leg_y = 0.5, 10.3
ax.text(leg_x, leg_y, 'Legend:', fontsize=8, fontweight='bold', color='#2c3e50')

# Activation arrow
arr_leg = FancyArrowPatch((leg_x, leg_y - 0.35), (leg_x + 0.8, leg_y - 0.35),
                           arrowstyle='->', mutation_scale=12, color=C_ACTIVATE, linewidth=1.5)
ax.add_patch(arr_leg)
ax.text(leg_x + 1.0, leg_y - 0.35, 'Activation', fontsize=7, va='center', color=C_ACTIVATE)

# Inhibition arrow
arr_leg2 = FancyArrowPatch((leg_x + 2.5, leg_y - 0.35), (leg_x + 3.3, leg_y - 0.35),
                            arrowstyle='-[', mutation_scale=8, color=C_INHIBIT, linewidth=1.5)
ax.add_patch(arr_leg2)
ax.text(leg_x + 3.5, leg_y - 0.35, 'Inhibition', fontsize=7, va='center', color=C_INHIBIT)

# Drug target
draw_node(leg_x + 5.5, leg_y - 0.35, 'Drug', color='#fadbd8', border=C_DRUG,
          fontsize=7, width=0.9, height=0.35)

# Dashed = indirect
ax.plot([leg_x + 6.5, leg_x + 7.3], [leg_y - 0.35, leg_y - 0.35],
        '--', color='#7f8c8d', linewidth=1.5)
ax.text(leg_x + 7.5, leg_y - 0.35, 'Indirect/multi-step', fontsize=7,
        va='center', color='#7f8c8d')

# =====================================================================
# PTCH1 COPY NUMBER ANNOTATION
# =====================================================================

ax.text(1.5, 6.55, 'Ptch1 CN: 1.0 (WT), 0.0 (MB)',
        ha='center', fontsize=7, color=C_HH_BORDER, fontstyle='italic',
        bbox=dict(boxstyle='round,pad=0.15', facecolor=C_HH, alpha=0.5))

# MYCN amplification annotation
ax.text(8.2, 7.55, 'Amp: 1x (WT), 2.8x (MB)',
        ha='center', fontsize=7, color=C_MYCN_BORDER, fontstyle='italic',
        bbox=dict(boxstyle='round,pad=0.15', facecolor=C_MYCN, alpha=0.5))

# =====================================================================
# SAVE
# =====================================================================

plt.tight_layout()
plt.savefig('simulations/fig_v42_architecture.png', dpi=250, bbox_inches='tight',
            facecolor='white')
plt.savefig('simulations/fig_v42_architecture.pdf', bbox_inches='tight',
            facecolor='white')
plt.close()
print("Saved: fig_v42_architecture.png and .pdf")
