from __future__ import annotations

import html
import json
import pandas as pd
import streamlit as st

# Form context for the main-area widgets (see the comment before the cost
# rail): pinned to streamlit==1.63.0, whose form internals these touch.
from streamlit.delta_generator_singletons import get_dg_singleton_instance
from streamlit.elements.lib.form_utils import FormData as _FormData

from economic_dna import (
    Scenario,
    find_crossover_years,
    simulate_dna_unit_costs,
    simulate_scenario,
    simulate_start_years,
)
from economic_dna.assumptions import load_assumptions
from economic_dna.visualization import (
    breakdown_chart,
    dna_unit_cost_chart,
    lifecycle_chart,
    palette_for,
    projection_chart,
)


st.set_page_config(
    page_title="DNA Storage Cost Explorer",
    layout="wide",
    # "auto" keeps the sidebar expanded on desktop and collapsed on phones,
    # where it renders as a full-screen overlay (see the mobile media query).
    initial_sidebar_state="auto",
)

_main_dg = get_dg_singleton_instance().main_dg
# The app temporarily attaches form data to the main dg so header/panel
# controls can join the sidebar form. A callback-triggered rerun can interrupt
# before the normal cleanup point, so each run must begin outside that form.
_main_dg._form_data = None

st.markdown(
    """
    <style>
    :root {
        --canvas: #eef2f0;
        --surface: #ffffff;
        --surface-subtle: #e6ece9;
        --sidebar-bg: #f8faf9;
        --ink: #18211f;
        --ink-soft: #34413d;
        --muted: #5f6e68;
        --line: #d9e0dc;
        --line-strong: #c2ccc7;
        --accent: #bd4b38;
        --accent-strong: #963829;
        --accent-soft: #faece8;
        --accent-glow: rgba(189, 75, 56, 0.32);
        --action-bg: #bd4b38;
        --action-bg-hover: #963829;
        --action-text: #ffffff;
        --teal: #1f716a;
        --gold: #b78322;
        --chip-bg: #f6ecd8;
        --chip-border: #e3cf9e;
        --chip-text: #8a6414;
        --download-text: #40504b;
        --scroll-thumb: #c2ccc7;
        --header-bg: rgba(238, 242, 240, 0.96);
        --shadow-sm: 0 1px 3px rgba(24, 33, 31, 0.07);
        --shadow-md: 0 3px 10px rgba(24, 33, 31, 0.14);
        --shadow-lg: 0 8px 22px rgba(24, 33, 31, 0.16);
        --ease: cubic-bezier(0.2, 0.7, 0.3, 1);
    }
    /* The dark palette applies when the Python-rendered theme marker carries
       theme-dark: the system preference decides the initial theme, and the
       sidebar toggle (or a shared URL) overrides it after. */
    html:has(.theme-marker.theme-dark) {
        --canvas: #0f1516;
        --surface: #171f20;
        --surface-subtle: #1e282a;
        --sidebar-bg: #121a1b;
        --ink: #e7edeb;
        --ink-soft: #c3cfcb;
        --muted: #93a19d;
        --line: #273234;
        --line-strong: #38464a;
        --accent: #d96a52;
        --accent-strong: #e58a74;
        --accent-soft: #3a231e;
        --accent-glow: rgba(217, 106, 82, 0.35);
        --action-bg: #c9533b;
        --action-bg-hover: #d96a52;
        --action-text: #ffffff;
        --teal: #45a694;
        --gold: #d4a63e;
        --chip-bg: #3a3018;
        --chip-border: #5a4a24;
        --chip-text: #e6c77e;
        --download-text: #aab9b4;
        --scroll-thumb: #38464a;
        --header-bg: rgba(15, 21, 22, 0.96);
        --shadow-sm: 0 1px 3px rgba(0, 0, 0, 0.45);
        --shadow-md: 0 3px 10px rgba(0, 0, 0, 0.55);
        --shadow-lg: 0 8px 22px rgba(0, 0, 0, 0.65);
        /* Streamlit theme tokens, so its own widgets follow the dark theme. */
        --background-color: #0f1516;
        --secondary-background-color: #171f20;
        --text-color: #e7edeb;
        --primary-color: #d96a52;
        --border-color: #273234;
    }

    /* Dark-mode fixes for Streamlit widgets whose emotion styles are compiled
       from the light theme config and ignore the CSS variables. */
    html:has(.theme-marker.theme-dark) .stButton button,
    html:has(.theme-marker.theme-dark) .stDownloadButton button {
        background: var(--surface);
        border-color: var(--line);
        color: var(--ink);
    }
    html:has(.theme-marker.theme-dark) .stButton button:hover,
    html:has(.theme-marker.theme-dark) .stDownloadButton button:hover {
        border-color: var(--accent);
        color: var(--accent-strong);
    }
    html:has(.theme-marker.theme-dark) [data-testid="stTooltipIcon"] {
        color: var(--accent) !important;
    }
    /* The expand-control icon is compiled from the light theme (60% dark ink,
       ignoring CSS vars) — invisible on the dark "Edit inputs" pill. */
    html:has(.theme-marker.theme-dark) [data-testid="stExpandSidebarButton"] [data-testid="stIconMaterial"] {
        color: var(--ink) !important;
        -webkit-text-fill-color: var(--ink) !important;
    }
    html:has(.theme-marker.theme-dark) [data-testid="stTooltipIcon"]:hover {
        color: var(--accent-strong) !important;
    }
    /* Streamlit's .icon class hardcodes a dark stroke; the "?" is drawn with
       stroke, not fill, so the stroke color itself must be overridden. */
    html:has(.theme-marker.theme-dark) [data-testid="stTooltipIcon"] svg {
        stroke: var(--accent) !important;
        stroke-width: 2.4;
    }
    /* The running indicator is drawn in dark ink and disappears on the dark
       header; brighten it to the dark theme's ink color. */
    html:has(.theme-marker.theme-dark) [data-testid="stStatusWidgetRunningIcon"] svg,
    html:has(.theme-marker.theme-dark) [data-testid="stStatusWidget"] span {
        color: var(--ink) !important;
    }
    html:has(.theme-marker.theme-dark) [data-testid="stStatusWidgetRunningIcon"] svg {
        filter: drop-shadow(0 0 4px rgba(231, 237, 235, 0.35));
    }
    html:has(.theme-marker.theme-dark) section[data-testid="stSidebar"] input[data-testid="stNumberInputField"],
    html:has(.theme-marker.theme-dark) section[data-testid="stSidebar"] [data-testid="stSelectbox"] input,
    html:has(.theme-marker.theme-dark) .st-key-cost-panel input[data-testid="stNumberInputField"],
    html:has(.theme-marker.theme-dark) .st-key-cost-panel input[data-testid="stTextInputField"] {
        color: var(--ink) !important;
        -webkit-text-fill-color: var(--ink) !important;
    }
    html:has(.theme-marker.theme-dark) section[data-testid="stSidebar"] [data-testid="stNumberInputStepUp"],
    html:has(.theme-marker.theme-dark) section[data-testid="stSidebar"] [data-testid="stNumberInputStepDown"],
    html:has(.theme-marker.theme-dark) .st-key-cost-panel [data-testid="stNumberInputStepUp"],
    html:has(.theme-marker.theme-dark) .st-key-cost-panel [data-testid="stNumberInputStepDown"] {
        color: var(--muted) !important;
    }
    html:has(.theme-marker.theme-dark) div:has(> [role="listbox"]) {
        background: var(--surface) !important;
        border: 1px solid var(--line-strong) !important;
    }
    html:has(.theme-marker.theme-dark) [role="listbox"] [role="option"] {
        background: transparent !important;
        color: var(--ink) !important;
    }
    html:has(.theme-marker.theme-dark) [role="listbox"] [role="option"]:hover,
    html:has(.theme-marker.theme-dark) [role="listbox"] [role="option"][aria-selected="true"] {
        background: var(--surface-subtle) !important;
    }

    html, body {
        font-family: "Aptos", "Segoe UI", Arial, sans-serif;
        letter-spacing: 0;
        /* Real phones inflate/boost font sizes on pages without this, which
           blows up fixed-size layouts and clips content. */
        -webkit-text-size-adjust: 100%;
        text-size-adjust: 100%;
    }
    .stApp {
        background: var(--canvas);
        color: var(--ink);
        font-family: "Aptos", "Segoe UI", Arial, sans-serif;
        transition: background-color 0.3s var(--ease), color 0.3s var(--ease);
    }
    [data-testid="stHeader"] {
        background: var(--header-bg);
        /* No backdrop-filter: at 0.96-alpha the blur is imperceptible, and the
           filter makes Chromium re-rasterize fixed elements layered above the
           header (the Calculate button) on every rerun — a visible flicker. */
        border-bottom: 1px solid var(--line);
        transition: background-color 0.3s var(--ease), border-color 0.3s var(--ease);
    }
    /* Theme toggle pinned to the page header, just left of Streamlit's toolbar
       menu (which occupies the rightmost 48px of the header). The running
       indicator sits immediately to the toggle's left instead of overlapping
       it. */
    .st-key-theme-toggle-anchor .stButton button {
        position: fixed;
        right: 3.9rem;
        top: 0.7rem;
        width: fit-content !important;
        z-index: 999990;
    }
    [data-testid="stStatusWidget"] {
        position: fixed !important;
        right: 9.9rem !important;
        top: 0.7rem !important;
        z-index: 999990;
    }
    .block-container {
        max-width: 1440px;
        /* Right gutter reserves space for the floating cost-assumption rail. */
        padding: 4.75rem 5.4rem 4rem 2.2rem;
    }

    h1, h2, h3, p { letter-spacing: 0; }
    h1 {
        color: var(--ink);
        font-size: 2.15rem;
        line-height: 1.14;
        margin: 0.15rem 0 0.35rem;
    }
    h2 { color: var(--ink); font-size: 1.25rem; }
    h3 { color: var(--ink); font-size: 1.05rem; }
    [data-testid="stCaptionContainer"] {
        color: var(--muted);
        font-size: 0.86rem;
        line-height: 1.55;
    }

    section[data-testid="stSidebar"] {
        width: 460px !important;
        background: var(--sidebar-bg);
        border-right: 1px solid var(--line);
        transition: background-color 0.3s var(--ease), border-color 0.3s var(--ease);
    }
    section[data-testid="stSidebar"] > div { width: 460px !important; }
    section[data-testid="stSidebar"] [data-testid="stSidebarContent"] {
        overflow-x: hidden;
        overflow-y: hidden;
        padding: 1.375rem 0.7rem 0 1.2rem;
    }
    section[data-testid="stSidebar"] [data-testid="stSidebarHeader"] {
        margin-bottom: 0;
    }
    section[data-testid="stSidebar"] [data-testid="stSidebarUserContent"] {
        padding-bottom: 0;
    }
    .sidebar-kicker, .page-kicker, .workspace-kicker, .export-label {
        color: var(--accent);
        font-size: 0.72rem;
        font-weight: 750;
        line-height: 1.2;
        letter-spacing: 0;
        text-transform: uppercase;
    }
    .sidebar-title {
        color: var(--ink);
        font-size: 1.35rem;
        font-weight: 720;
        line-height: 1.25;
        margin: 0.22rem 0 0.35rem;
    }
    .sidebar-copy {
        color: var(--muted);
        font-size: 0.86rem;
        line-height: 1.5;
        margin: 0 0 0.9rem;
    }
    .sidebar-section {
        border-top: 1px solid var(--line);
        color: var(--ink);
        font-size: 0.77rem;
        font-weight: 720;
        letter-spacing: 0;
        margin: 1.1rem 0 0.65rem;
        padding-top: 0.9rem;
        text-transform: uppercase;
    }
    section[data-testid="stSidebar"] label p,
    .st-key-cost-panel label p {
        color: var(--ink-soft);
        font-size: 0.82rem;
        font-weight: 600;
    }
    .log-scale-label {
        color: var(--ink-soft);
        font-size: 0.82rem;
        font-weight: 600;
        margin: 0 0 0.3rem;
    }
    section[data-testid="stSidebar"] [data-testid="stNumberInputContainer"],
    section[data-testid="stSidebar"] [data-testid="stSelectbox"] div:has(> input),
    .st-key-cost-panel [data-testid="stNumberInputContainer"] {
        background: var(--surface);
        border-color: var(--line-strong);
        border-radius: 6px;
    }
    /* Text inputs render borderless in 1.63 — the visible border goes on the
       INNER root element so it contains only the text box, not the label
       (both themes via the vars). */
    .st-key-cost-panel [data-testid="stTextInputRootElement"] {
        background: var(--surface);
        border: 1px solid var(--line-strong);
        border-radius: 6px;
    }
    section[data-testid="stSidebar"] [data-testid="stNumberInputContainer"]:focus-within,
    section[data-testid="stSidebar"] [data-testid="stSelectbox"] div:has(> input):focus-within,
    .st-key-cost-panel [data-testid="stNumberInputContainer"]:focus-within {
        border-color: var(--accent);
        box-shadow: 0 0 0 1px var(--accent);
    }
    /* The focus ring goes on the text BOX, not the wrapper (which also
       contains the label — a ring there would outline "Display name" too). */
    .st-key-cost-panel [data-testid="stTextInputRootElement"]:focus-within {
        border-color: var(--accent);
        box-shadow: 0 0 0 1px var(--accent);
    }
    section[data-testid="stSidebar"] [data-testid="stVerticalBlockBorderWrapper"] {
        background: var(--surface);
        border-color: var(--line);
        border-radius: 6px;
    }

    .stButton button, .stDownloadButton button {
        border-radius: 8px;
        font-weight: 650;
        min-height: 2.35rem;
        transition: border-color 120ms ease, background-color 120ms ease, color 120ms ease,
            transform 90ms ease, box-shadow 120ms ease;
    }
    .stButton button:hover, .stDownloadButton button:hover {
        border-color: var(--accent);
        color: var(--accent-strong);
    }
    .stButton button:active, .stDownloadButton button:active {
        transform: translateY(1px) scale(0.985);
    }
    [data-testid="stFormSubmitButton"] button {
        min-height: 2.75rem;
        background: var(--accent);
        border-color: var(--accent);
        color: #ffffff;
        font-weight: 700;
    }
    [data-testid="stFormSubmitButton"] button:hover {
        background: var(--accent-strong);
        border-color: var(--accent-strong);
        color: #ffffff;
    }
    /* Two-column sidebar: the left column holds all inputs and scrolls on its
       own; the right column holds only the Calculate button. The sidebar's
       own scroll container is disabled (position: sticky is unreliable inside
       it in this DOM: the absolutely-positioned flex app container breaks it),
       so instead the input column becomes the scroller and the button rail
       sits beside it as a fixed sibling that never moves. The in-flow sidebar
       header occupies the top 82px (22px content padding + 60px header in
       Streamlit 1.63), which is reserved so nothing overlaps the collapse
       control. */
    section[data-testid="stSidebar"] [data-testid="stHorizontalBlock"]:has(.st-key-scenario-action-rail) {
        align-items: stretch;
        gap: 0.7rem !important;
    }
    /* Note the child combinator: a descendant selector would also match the
       first column of every inner st.columns() row inside the input column,
       stretching each widget row to the full visible height. */
    section[data-testid="stSidebar"] [data-testid="stHorizontalBlock"]:has(.st-key-scenario-action-rail)
    > [data-testid="stColumn"]:first-child {
        height: calc(100vh - 82px);
        overflow-y: auto;
        overscroll-behavior-y: contain;
        padding: 0 0 1.5rem;
        scrollbar-gutter: stable;
    }
    section[data-testid="stSidebar"] [data-testid="stHorizontalBlock"]:has(.st-key-scenario-action-rail)
    > [data-testid="stColumn"]:first-child::-webkit-scrollbar {
        width: 8px;
    }
    section[data-testid="stSidebar"] [data-testid="stHorizontalBlock"]:has(.st-key-scenario-action-rail)
    > [data-testid="stColumn"]:first-child::-webkit-scrollbar-thumb {
        background: var(--scroll-thumb);
        border-radius: 4px;
    }
    section[data-testid="stSidebar"] .st-key-scenario-action-rail {
        flex: 0 0 auto;
        height: calc(100vh - 82px);
        min-height: calc(100vh - 82px);
        overflow: visible;
        position: relative;
        z-index: 1;
    }
    /* Percentage heights cannot resolve here (the rail sits in a flex-content-
       sized column), so the button fills the rail with absolute positioning
       instead of a height: 100% chain. */
    section[data-testid="stSidebar"] .st-key-scenario-action-rail
    [data-testid="stElementContainer"] {
        position: absolute;
        inset: 0;
    }
    @keyframes pending-pulse {
        0%, 100% { box-shadow: var(--shadow-md), inset 0 1px 0 rgba(255, 255, 255, 0.22),
            0 0 0 0 var(--accent-glow); }
        55% { box-shadow: var(--shadow-md), inset 0 1px 0 rgba(255, 255, 255, 0.22),
            0 0 0 6px transparent; }
    }
    section[data-testid="stSidebar"] .st-key-scenario-action-rail [data-testid="stFormSubmitButton"] button {
        align-items: center;
        background: linear-gradient(168deg, var(--action-bg) 0%, var(--action-bg-hover) 130%);
        border-color: var(--action-bg);
        border-radius: 8px;
        box-shadow: var(--shadow-md);
        color: var(--action-text);
        display: flex;
        flex-direction: column;
        inset: 0;
        justify-content: space-between;
        padding: 1rem 0.45rem;
        position: absolute;
        transition: filter 140ms ease, transform 90ms ease, box-shadow 140ms ease;
    }
    section[data-testid="stSidebar"] .st-key-scenario-action-rail [data-testid="stFormSubmitButton"] button:hover {
        filter: brightness(1.09);
        transform: translateY(-1px);
        box-shadow: var(--shadow-lg);
    }
    section[data-testid="stSidebar"] .st-key-scenario-action-rail [data-testid="stFormSubmitButton"] button:active {
        filter: brightness(0.94);
        transform: translateY(1px) scale(0.992);
        box-shadow: var(--shadow-sm);
    }
    .pending-marker { display: none; }
    [data-testid="stMarkdownContainer"]:has(.theme-marker) { display: none; }
    .st-key-theme_auto_dark,
    [data-testid="stElementContainer"]:has(.st-key-theme_auto_dark) {
        display: none;
    }
    /* The 1x1 utility iframes (theme sync, click-away) would otherwise render
       as tiny bordered dashes in the main flow. Hidden, but their scripts
       still load and run. */
    [data-testid="stIFrame"] {
        display: none;
    }
    section[data-testid="stSidebar"] [data-testid="stSidebarCollapseButton"],
    [data-testid="stExpandSidebarButton"],
    [data-testid="stSidebarCollapsedControl"] {
        pointer-events: auto !important;
        z-index: 5 !important;
    }
    /* When the sidebar is collapsed its content slides off-canvas but still
       overlaps the main area (the sidebar section keeps a high z-index), so
       the pinned button would paint over the main content and the header
       would intercept clicks aimed at the expand control. Hide the rail and
       let pointer events pass through to the main area. */
    section[data-testid="stSidebar"][aria-expanded="false"] .st-key-scenario-action-rail {
        visibility: hidden;
    }
    section[data-testid="stSidebar"][aria-expanded="false"] [data-testid="stSidebarContent"] {
        pointer-events: none;
    }
    section[data-testid="stSidebar"] .st-key-scenario-action-rail
    [data-testid="stFormSubmitButton"] button::before,
    section[data-testid="stSidebar"] .st-key-scenario-action-rail
    [data-testid="stFormSubmitButton"] button::after {
        content: "\\2192";
        display: block;
        flex: 0 0 auto;
        font-family: "Segoe UI Symbol", sans-serif;
        font-size: 1.55rem;
        line-height: 1;
    }
    section[data-testid="stSidebar"] .st-key-scenario-action-rail [data-testid="stFormSubmitButton"] button p {
        font-size: 1.3rem;
        font-weight: 700;
        line-height: 1.2;
        white-space: nowrap;
        writing-mode: vertical-rl;
        transform: rotate(180deg);
    }
    .stDownloadButton button {
        background: var(--surface);
        border-color: var(--line);
        color: var(--download-text);
    }
    [class*="st-key-download_"] .stButton button,
    [class*="st-key-download_"] button {
        background: var(--surface);
        border-color: var(--line);
        color: var(--download-text);
    }
    .chart-loading {
        align-items: center;
        background:
            linear-gradient(90deg, transparent 0%, color-mix(in srgb, var(--surface-subtle) 70%, transparent) 48%, transparent 96%),
            var(--surface);
        background-size: 220% 100%, 100% 100%;
        border: 1px solid var(--line);
        border-radius: 8px;
        box-sizing: border-box;
        display: flex;
        justify-content: center;
        margin: 0.35rem 0 0.85rem;
        min-height: 240px;
        overflow: hidden;
        position: relative;
        width: 100%;
        animation: chart-loading-sheen 1.25s var(--ease) infinite;
    }
    .chart-loading::before {
        color: var(--muted);
        content: "Loading chart";
        font-size: 0.82rem;
        font-weight: 650;
    }
    @keyframes chart-loading-sheen {
        0% { background-position: 150% 0, 0 0; }
        100% { background-position: -70% 0, 0 0; }
    }

    /* On desktop the sidebar's Calculate is hidden: the wide centered header
       button is the single always-visible Calculate (phones keep the bar,
       see the 640px media query). */
    section[data-testid="stSidebar"] .st-key-scenario-action-rail {
        display: none;
    }
    /* Header Calculate: pinned to the top bar, centered in the space between
       the left controls and the theme toggle / status widget on the right,
       with downward arrows flanking the label. It pulses while changes are
       pending. */
    .st-key-calculate-anchor [data-testid="stFormSubmitButton"] button {
        background: linear-gradient(168deg, var(--action-bg) 0%, var(--action-bg-hover) 130%);
        border: none;
        border-radius: 12px;
        box-shadow: var(--shadow-md), inset 0 1px 0 rgba(255, 255, 255, 0.22);
        color: var(--action-text);
        font-size: 1.05rem;
        font-weight: 700;
        gap: 0.55rem;
        left: 50%;
        letter-spacing: 0.03em;
        min-height: 2.6rem;
        position: fixed;
        top: 0.7rem;
        transform: translateX(-50%) translateZ(0);
        width: min(44rem, calc(100vw - 26rem)) !important;
        z-index: 999990;
        /* Own compositing layer: reruns repaint the content below, and a
           promoted layer keeps the button from blinking along with it. */
        will-change: transform;
    }
    .st-key-calculate-anchor [data-testid="stFormSubmitButton"] button::before,
    .st-key-calculate-anchor [data-testid="stFormSubmitButton"] button::after {
        content: "↓";
        font-size: 1.05rem;
        font-weight: 700;
        line-height: 1;
    }
    .st-key-calculate-anchor [data-testid="stFormSubmitButton"] button:hover {
        box-shadow: var(--shadow-lg), inset 0 1px 0 rgba(255, 255, 255, 0.22);
        filter: brightness(1.07);
        transform: translateX(-50%) translateZ(0) translateY(-1px);
    }
    .st-key-calculate-anchor [data-testid="stFormSubmitButton"] button:active {
        transform: translateX(-50%) translateZ(0) translateY(1px) scale(0.99);
    }
    /* During a rerun Streamlit disables buttons; keep the header button
       visually identical so the disabled state never reads as a blink. */
    .st-key-calculate-anchor [data-testid="stFormSubmitButton"] button[disabled],
    .st-key-calculate-anchor [data-testid="stFormSubmitButton"] button:disabled {
        background: linear-gradient(168deg, var(--action-bg) 0%, var(--action-bg-hover) 130%) !important;
        color: var(--action-text) !important;
        opacity: 1 !important;
    }
    html:has(.theme-marker.theme-dark) .st-key-calculate-anchor [data-testid="stFormSubmitButton"] button {
        background: linear-gradient(168deg, var(--action-bg) 0%, var(--action-bg-hover) 130%);
        color: var(--action-text);
    }
    /* Radio tab text: the option's inner text div is compiled from the light
       theme (dark ink) and ignores the label color — inherit it instead. */
    html:has(.theme-marker.theme-dark) .st-key-cost-rail [role="radiogroup"] label div {
        color: inherit !important;
    }
    /* The custom-name text input can paint a light-theme (or autofill) white
       box inside the dark panel — force the dark surface on the input itself,
       including Chrome's autofill overlay. */
    html:has(.theme-marker.theme-dark) .st-key-cost-panel [data-testid="stTextInput"] input {
        background: var(--surface) !important;
        -webkit-box-shadow: inset 0 0 0 1000px var(--surface) !important;
    }
    html:has(.theme-marker.theme-dark) .st-key-cost-panel [data-testid="stTextInput"] input:-webkit-autofill {
        -webkit-box-shadow: inset 0 0 0 1000px var(--surface) !important;
        -webkit-text-fill-color: var(--ink) !important;
    }
    /* The panel's Calculate bar (phones only — the header button covers
       desktop) sits sticky at the bottom of the scrolling panel sheet. */
    .st-key-calculate_panel {
        display: none;
    }
    html:has(.theme-marker.theme-dark) .st-key-calculate_panel [data-testid="stFormSubmitButton"] button {
        background: var(--action-bg);
        border-color: var(--action-bg);
        color: var(--action-text);
    }

    /* Cost-assumption rail: a compact model switcher pinned to the edge.
       The panel stays mounted and animates with opacity/transform so both
       opening and closing feel smooth while form edits stay intact. */
    .st-key-cost-rail {
        background: color-mix(in srgb, var(--surface) 94%, transparent);
        border: 1px solid var(--line);
        border-radius: 16px;
        box-shadow: var(--shadow-md);
        padding: 0.38rem;
        position: fixed;
        right: 0.75rem;
        top: 5.05rem;
        z-index: 999900;
        width: 8.5rem;
    }
    .st-key-cost-rail [data-testid="stVerticalBlock"] {
        gap: 0.28rem !important;
    }
    .st-key-cost-rail [data-testid="stElementContainer"] {
        width: 100% !important;
    }
    .st-key-cost-rail [role="radiogroup"] {
        flex-direction: column;
        gap: 0.28rem !important;
        width: 100%;
    }
    .st-key-cost-rail [role="radiogroup"] label {
        align-items: center;
        background: transparent;
        border: 1px solid transparent;
        border-radius: 11px;
        box-shadow: none;
        color: var(--muted);
        cursor: pointer;
        display: flex;
        box-sizing: border-box;
        font-size: 0.7rem;
        font-weight: 700;
        gap: 0.4rem;
        justify-content: flex-start;
        letter-spacing: 0;
        margin: 0 !important;
        min-height: 2.7rem;
        overflow: hidden;
        padding: 0.28rem 0.36rem;
        position: relative;
        text-transform: none;
        transform: none;
        transition: border-color 160ms var(--ease), background-color 160ms var(--ease),
            color 160ms var(--ease), transform 120ms var(--ease), box-shadow 160ms var(--ease);
        user-select: none;
        white-space: nowrap;
        width: 100%;
        writing-mode: horizontal-tb;
    }
    /* The ✕ option is the radio's "closed" state and must remain in the group,
       but it is completely invisible: the panel closes by clicking the open
       tab again or clicking away, so no visible close button is needed. The
       synthesized closeLabel.click() in the click-away script still works on
       a display:none label. */
    .st-key-cost-rail [role="radiogroup"] label:has(input[value="0"]) {
        display: none;
    }
    .st-key-cost-rail [role="radiogroup"] label:has(input:focus-visible) {
        outline: 2px solid var(--accent);
        outline-offset: 2px;
    }
    /* Compact badges make each storage model scannable without relying only
       on the text label. */
    .st-key-cost-rail [role="radiogroup"] label::before {
        align-items: center;
        background: var(--surface-subtle);
        border: 1px solid var(--line);
        border-radius: 999px;
        color: var(--ink-soft);
        content: "";
        display: flex;
        flex: 0 0 1.45rem;
        font-size: 0.64rem;
        font-weight: 800;
        height: 1.45rem;
        justify-content: center;
        line-height: 1;
        transition: background-color 160ms var(--ease), border-color 160ms var(--ease),
            color 160ms var(--ease), transform 160ms var(--ease);
        width: 1.45rem;
    }
    .st-key-cost-rail [role="radiogroup"] label:has(input[value="1"])::before { content: "DNA"; }
    .st-key-cost-rail [role="radiogroup"] label:has(input[value="2"])::before { content: "S3"; }
    .st-key-cost-rail [role="radiogroup"] label:has(input[value="3"])::before { content: "AZ"; }
    .st-key-cost-rail [role="radiogroup"] label:has(input[value="4"])::before { content: "LTO"; }
    .st-key-cost-rail [role="radiogroup"] label:has(input[value="5"])::before { content: "+"; }
    .st-key-cost-rail [role="radiogroup"] label::after {
        background: var(--accent);
        border-radius: 999px;
        content: "";
        height: 0.36rem;
        margin-left: auto;
        opacity: 0;
        transform: scale(0.4);
        transition: opacity 160ms var(--ease), transform 160ms var(--ease);
        width: 0.36rem;
    }
    .st-key-cost-rail [role="radiogroup"] label:hover {
        background: var(--surface-subtle);
        border-color: var(--accent);
        color: var(--accent-strong);
        transform: translateX(-2px);
    }
    .st-key-cost-rail [role="radiogroup"] label:active {
        transform: translateX(-1px) scale(0.985);
    }
    /* The open tab flattens against the panel: same surface, no left corner,
       accent ink and an accent bar mark the active model. */
    body:has(.st-key-cost_model_radio input[value="1"]:checked) .st-key-cost-rail label:has(input[value="1"]),
    body:has(.st-key-cost_model_radio input[value="2"]:checked) .st-key-cost-rail label:has(input[value="2"]),
    body:has(.st-key-cost_model_radio input[value="3"]:checked) .st-key-cost-rail label:has(input[value="3"]),
    body:has(.st-key-cost_model_radio input[value="4"]:checked) .st-key-cost-rail label:has(input[value="4"]),
    body:has(.st-key-cost_model_radio input[value="5"]:checked) .st-key-cost-rail label:has(input[value="5"]) {
        background: var(--accent-soft);
        border-color: color-mix(in srgb, var(--accent) 45%, var(--line));
        color: var(--accent-strong);
        box-shadow: inset 3px 0 0 var(--accent);
        transform: translateX(-3px);
    }
    body:has(.st-key-cost_model_radio input[value="1"]:checked) .st-key-cost-rail label:has(input[value="1"])::before,
    body:has(.st-key-cost_model_radio input[value="2"]:checked) .st-key-cost-rail label:has(input[value="2"])::before,
    body:has(.st-key-cost_model_radio input[value="3"]:checked) .st-key-cost-rail label:has(input[value="3"])::before,
    body:has(.st-key-cost_model_radio input[value="4"]:checked) .st-key-cost-rail label:has(input[value="4"])::before,
    body:has(.st-key-cost_model_radio input[value="5"]:checked) .st-key-cost-rail label:has(input[value="5"])::before {
        background: var(--accent);
        border-color: var(--accent);
        color: #ffffff;
        transform: scale(1.04);
    }
    body:has(.st-key-cost_model_radio input[value="1"]:checked) .st-key-cost-rail label:has(input[value="1"])::after,
    body:has(.st-key-cost_model_radio input[value="2"]:checked) .st-key-cost-rail label:has(input[value="2"])::after,
    body:has(.st-key-cost_model_radio input[value="3"]:checked) .st-key-cost-rail label:has(input[value="3"])::after,
    body:has(.st-key-cost_model_radio input[value="4"]:checked) .st-key-cost-rail label:has(input[value="4"])::after,
    body:has(.st-key-cost_model_radio input[value="5"]:checked) .st-key-cost-rail label:has(input[value="5"])::after {
        opacity: 1;
        transform: scale(1);
    }
    .st-key-cost-panel {
        background: var(--surface);
        border: 1px solid var(--line);
        border-radius: 16px;
        box-shadow: var(--shadow-lg);
        /* Hidden by default: during initial hydration the radio has no
           selection yet, and the :checked-based hide rule cannot match —
           a default-visible panel flashes for those first frames. */
        display: flex;
        flex-direction: column;
        max-height: calc(100vh - 7rem);
        opacity: 0;
        overflow-y: auto;
        pointer-events: none;
        padding: 1.1rem 1.15rem 1.2rem;
        position: fixed;
        right: 9.85rem;
        top: 5.05rem;
        transform: translateX(1.2rem) scale(0.985);
        transition: opacity 220ms var(--ease), transform 260ms var(--ease),
            visibility 0s linear 260ms, border-color 180ms var(--ease);
        visibility: hidden;
        width: min(27rem, calc(100vw - 11.5rem));
        z-index: 999900;
    }
    body:has(.st-key-cost_model_radio input[value="1"]:checked) .st-key-cost-panel,
    body:has(.st-key-cost_model_radio input[value="2"]:checked) .st-key-cost-panel,
    body:has(.st-key-cost_model_radio input[value="3"]:checked) .st-key-cost-panel,
    body:has(.st-key-cost_model_radio input[value="4"]:checked) .st-key-cost-panel,
    body:has(.st-key-cost_model_radio input[value="5"]:checked) .st-key-cost-panel {
        opacity: 1;
        pointer-events: auto;
        transform: none;
        transition-delay: 0s;
        visibility: visible;
    }
    .st-key-cost-panel [data-testid="stLayoutWrapper"]:has(.st-key-cost_model_dna),
    .st-key-cost-panel [data-testid="stLayoutWrapper"]:has(.st-key-cost_model_amazon),
    .st-key-cost-panel [data-testid="stLayoutWrapper"]:has(.st-key-cost_model_azure),
    .st-key-cost-panel [data-testid="stLayoutWrapper"]:has(.st-key-cost_model_tape),
    .st-key-cost-panel [data-testid="stLayoutWrapper"]:has(.st-key-cost_model_custom) {
        display: none;
    }
    body:has(.st-key-cost_model_radio input[value="1"]:checked) .st-key-cost-panel [data-testid="stLayoutWrapper"]:has(.st-key-cost_model_dna),
    body:has(.st-key-cost_model_radio input[value="2"]:checked) .st-key-cost-panel [data-testid="stLayoutWrapper"]:has(.st-key-cost_model_amazon),
    body:has(.st-key-cost_model_radio input[value="3"]:checked) .st-key-cost-panel [data-testid="stLayoutWrapper"]:has(.st-key-cost_model_azure),
    body:has(.st-key-cost_model_radio input[value="4"]:checked) .st-key-cost-panel [data-testid="stLayoutWrapper"]:has(.st-key-cost_model_tape),
    body:has(.st-key-cost_model_radio input[value="5"]:checked) .st-key-cost-panel [data-testid="stLayoutWrapper"]:has(.st-key-cost_model_custom) {
        display: block;
        animation: cost-model-in 210ms var(--ease) both;
    }
    @keyframes cost-model-in {
        from { opacity: 0; transform: translateY(6px); }
        to { opacity: 1; transform: none; }
    }
    /* The panel title follows the checked tab; the default title shows when
       the panel is closed. */
    .st-key-cost-panel .cost-panel-title { display: none; }
    body:has(.st-key-cost_model_radio input[value="0"]:checked) .st-key-cost-panel .cost-title-none { display: block; }
    body:has(.st-key-cost_model_radio input[value="1"]:checked) .st-key-cost-panel .cost-title-dna { display: block; }
    body:has(.st-key-cost_model_radio input[value="2"]:checked) .st-key-cost-panel .cost-title-amazon { display: block; }
    body:has(.st-key-cost_model_radio input[value="3"]:checked) .st-key-cost-panel .cost-title-azure { display: block; }
    body:has(.st-key-cost_model_radio input[value="4"]:checked) .st-key-cost-panel .cost-title-tape { display: block; }
    body:has(.st-key-cost_model_radio input[value="5"]:checked) .st-key-cost-panel .cost-title-custom { display: block; }

    .st-key-cost-panel .cost-panel-title {
        color: var(--ink);
        font-size: 1.2rem;
        font-weight: 720;
        line-height: 1.3;
        margin: 0.15rem 2.4rem 0.35rem 0;
    }
    .st-key-cost-panel .cost-panel-deck {
        color: var(--muted);
        font-size: 0.84rem;
        line-height: 1.5;
        margin: 0;
    }
    .st-key-cost-panel .cost-panel-header {
        background: linear-gradient(180deg, var(--surface-subtle), transparent 115%);
        border-bottom: 1px solid var(--line);
        margin: -1.1rem -1.15rem 1rem;
        padding: 1.05rem 1.15rem 0.95rem;
        position: sticky;
        top: -1.1rem;
        z-index: 1;
    }
    .st-key-cost-panel .cost-panel-summary {
        display: none;
        flex-wrap: wrap;
        gap: 0.35rem;
        margin-top: 0.7rem;
    }
    body:has(.st-key-cost_model_radio input[value="1"]:checked) .st-key-cost-panel .cost-summary-dna,
    body:has(.st-key-cost_model_radio input[value="2"]:checked) .st-key-cost-panel .cost-summary-amazon,
    body:has(.st-key-cost_model_radio input[value="3"]:checked) .st-key-cost-panel .cost-summary-azure,
    body:has(.st-key-cost_model_radio input[value="4"]:checked) .st-key-cost-panel .cost-summary-tape,
    body:has(.st-key-cost_model_radio input[value="5"]:checked) .st-key-cost-panel .cost-summary-custom {
        display: flex;
    }
    .st-key-cost-panel .cost-panel-summary span {
        background: var(--surface);
        border: 1px solid var(--line);
        border-radius: 999px;
        color: var(--ink-soft);
        font-size: 0.72rem;
        font-weight: 650;
        line-height: 1.2;
        padding: 0.28rem 0.5rem;
    }
    .st-key-cost-panel [data-testid="stCaptionContainer"] {
        background: var(--surface-subtle);
        border: 1px solid var(--line);
        border-radius: 7px;
        color: var(--teal);
        font-size: 0.72rem;
        font-weight: 750;
        line-height: 1.2;
        margin: 0.75rem 0 0.35rem;
        padding: 0.34rem 0.5rem;
        text-transform: uppercase;
    }
    .st-key-cost-panel [data-testid="stNumberInput"],
    .st-key-cost-panel [data-testid="stTextInput"] {
        margin-bottom: 0.18rem;
    }
    /* Reset buttons synchronize Streamlit's widget state without updating the
       committed scenario or refreshing charts. */
    .st-key-sidebar-reset-btn .stButton button {
        background: transparent;
        border: 1px solid var(--line);
        border-radius: 8px;
        color: var(--muted);
        cursor: pointer;
        display: block;
        font-size: 0.9rem;
        font-weight: 650;
        min-height: 2.35rem;
        padding: 0 1rem;
        width: 100%;
    }
    .st-key-sidebar-reset-btn .stButton button:hover {
        background: var(--surface-subtle);
        border-color: var(--accent);
        color: var(--accent-strong);
    }
    .st-key-global-reset-btn .stButton button {
        background: var(--surface);
        border: 1px solid var(--line-strong);
        border-radius: 8px;
        box-shadow: var(--shadow-sm);
        color: var(--download-text);
        cursor: pointer;
        font-size: 0.88rem;
        font-weight: 650;
        left: calc(50% + min(22rem, calc((100vw - 26rem) / 2)) + 0.75rem);
        min-height: 2.35rem;
        padding: 0 1rem;
        position: fixed;
        top: 0.82rem;
        white-space: nowrap;
        width: max-content !important;
        z-index: 999990;
    }
    .st-key-global-reset-btn .stButton button:hover {
        background: var(--surface-subtle);
        border-color: var(--accent);
        color: var(--accent-strong);
    }
    /* Per-model resets sit inside the active cost panel and restore only that
       model's editable assumptions. */
    .st-key-cost-panel [class*="st-key-model-reset-"] [data-testid="stFormSubmitButton"] {
        display: flex;
        justify-content: flex-end;
    }
    .st-key-cost-panel [class*="st-key-model-reset-"] [data-testid="stFormSubmitButton"] button {
        background: transparent;
        border: 1px solid var(--line);
        border-radius: 999px;
        color: var(--muted);
        cursor: pointer;
        display: block;
        font-size: 0.85rem;
        font-weight: 650;
        margin: 0 0 0.7rem auto;
        min-height: 2.4rem;
        padding: 0 1.1rem;
        width: auto;
    }
    .st-key-cost-panel [class*="st-key-model-reset-"] [data-testid="stFormSubmitButton"] button:hover {
        background: var(--surface-subtle);
        border-color: var(--accent);
        color: var(--accent-strong);
    }

    .page-kicker { margin-bottom: 0.35rem; }
    .page-deck {
        color: var(--muted);
        font-size: 1rem;
        line-height: 1.55;
        margin: 0 0 1.15rem;
        max-width: 760px;
    }
    .model-strip {
        align-items: center;
        background: var(--surface);
        border: 1px solid var(--line);
        border-left: 3px solid var(--accent);
        border-radius: 12px;
        color: var(--muted);
        display: flex;
        flex-wrap: wrap;
        font-size: 0.8rem;
        gap: 0.35rem 0;
        margin: 0 0 1.25rem;
        min-height: 42px;
        padding: 0.65rem 0.85rem;
        transition: background-color 0.3s var(--ease), border-color 0.3s var(--ease);
    }
    .model-strip strong { color: var(--ink); font-weight: 700; }
    .model-strip a {
        color: var(--accent-strong);
        font-weight: 700;
        text-decoration: none;
    }
    .model-strip a:hover { text-decoration: underline; }
    .model-item {
        border-right: 1px solid var(--line);
        margin-right: 0.75rem;
        padding-right: 0.75rem;
    }
    .model-item:last-child { border-right: 0; margin-right: 0; padding-right: 0; }

    .scenario-bar {
        align-items: center;
        background: var(--surface-subtle);
        border: 1px solid var(--line);
        border-radius: 12px;
        color: var(--muted);
        display: flex;
        flex-wrap: wrap;
        font-size: 0.82rem;
        gap: 0.45rem 1.15rem;
        margin: 0 0 1rem;
        min-height: 44px;
        padding: 0.65rem 0.85rem;
        transition: background-color 0.3s var(--ease), border-color 0.3s var(--ease);
    }
    .scenario-bar strong { color: var(--ink); }
    .scenario-label { color: var(--teal); font-weight: 750; text-transform: uppercase; }
    @keyframes chip-in {
        from { opacity: 0; transform: translateX(-5px); }
        to { opacity: 1; transform: none; }
    }
    @keyframes fade-up {
        from { opacity: 0; transform: translateY(7px); }
        to { opacity: 1; transform: none; }
    }
    [data-testid="stMetric"] {
        animation: fade-up 0.45s var(--ease) both;
        background: var(--surface);
        border: 1px solid var(--line);
        border-radius: 12px;
        box-shadow: var(--shadow-sm);
        min-height: 114px;
        overflow: hidden;
        padding: 0.95rem 1rem;
        position: relative;
        transition: transform 160ms var(--ease), box-shadow 160ms var(--ease),
            background-color 0.3s var(--ease), border-color 0.3s var(--ease);
    }
    [data-testid="stMetric"]:hover {
        box-shadow: var(--shadow-md);
        transform: translateY(-2px);
    }
    [data-testid="stMetric"]::before {
        background: linear-gradient(90deg, var(--accent) 0%, var(--gold) 130%);
        content: "";
        height: 3px;
        left: 0;
        position: absolute;
        right: 0;
        top: 0;
    }
    [data-testid="stMetricLabel"] p {
        color: var(--muted);
        font-size: 0.79rem;
        font-weight: 650;
        line-height: 1.35;
        white-space: normal;
    }
    [data-testid="stMetricValue"] {
        color: var(--ink);
        font-size: 1.48rem;
        font-weight: 700;
        line-height: 1.2;
    }
    button:focus-visible, input:focus-visible {
        outline: 2px solid var(--accent);
        outline-offset: 2px;
    }
    [data-testid="stMetricDelta"] { color: var(--muted); }

    .workspace-kicker { margin: 1.55rem 0 0.35rem; }
    div[data-testid="stTabs"] div[role="tablist"] {
        background: var(--surface);
        border: 1px solid var(--line);
        border-radius: 12px;
        box-shadow: var(--shadow-sm);
        gap: 0.3rem;
        margin: 0 0 1.2rem;
        max-width: 100%;
        overflow-x: auto;
        padding: 0.3rem;
        width: fit-content;
    }
    div[data-testid="stTabs"] div[data-testid="stTab"][role="tab"] {
        background: transparent;
        border: 0;
        border-radius: 9px;
        color: var(--muted);
        flex: 0 0 auto;
        font-weight: 650;
        min-height: 48px;
        padding: 0.7rem 1.35rem;
        transition: background-color 140ms ease, color 140ms ease;
    }
    div[data-testid="stTabs"] div[data-testid="stTab"][role="tab"] p {
        font-size: 1rem;
    }
    div[data-testid="stTabs"] div[data-testid="stTab"][role="tab"]:hover {
        background: var(--surface-subtle);
        color: var(--ink);
    }
    div[data-testid="stTabs"] div[data-testid="stTab"][role="tab"][aria-selected="true"] {
        background: var(--action-bg);
        color: var(--action-text);
        font-weight: 700;
    }
    div[data-testid="stTabs"] div[data-testid="stTab"][role="tab"][aria-selected="true"] p {
        color: var(--action-text);
    }
    div[data-testid="stTabs"] .react-aria-SelectionIndicator,
    div[data-testid="stTabs"] div[data-baseweb="tab-highlight"],
    div[data-testid="stTabs"] div[data-baseweb="tab-border"] { display: none; }

    .tab-intro { margin: 0.15rem 0 1rem; max-width: 900px; }
    .tab-intro h2 {
        font-size: 1.18rem;
        line-height: 1.3;
        margin: 0 0 0.25rem;
    }
    .tab-intro p {
        color: var(--muted);
        font-size: 0.88rem;
        line-height: 1.55;
        margin: 0;
    }
    .chart-divider {
        border-top: 1px solid var(--line);
        margin: 1.65rem 0 1rem;
        padding-top: 1.15rem;
    }
    .chart-divider h3 { font-size: 1rem; margin: 0 0 0.2rem; }
    .chart-divider p { color: var(--muted); font-size: 0.84rem; margin: 0; }
    div[data-testid="stPlotlyChart"] {
        animation: fade-up 0.5s var(--ease) both;
        background: var(--surface);
        border: 1px solid var(--line);
        border-radius: 12px;
        box-shadow: var(--shadow-sm);
        overflow: hidden;
        padding: 0.2rem;
        transition: background-color 0.3s var(--ease), border-color 0.3s var(--ease),
            box-shadow 160ms var(--ease);
    }
    div[data-testid="stPlotlyChart"]:hover {
        box-shadow: var(--shadow-md);
    }
    .contract-card {
        background: var(--surface);
        border: 1px solid var(--line);
        border-radius: 12px;
        box-shadow: var(--shadow-sm);
        overflow: hidden;
        transition: background-color 0.3s var(--ease), border-color 0.3s var(--ease);
    }
    .contract-table {
        border-collapse: collapse;
        font-size: 0.88rem;
        width: 100%;
    }
    .contract-table th, .contract-table td {
        border-bottom: 1px solid var(--line);
        padding: 0.55rem 0.8rem;
        text-align: left;
        vertical-align: top;
    }
    .contract-table th {
        color: var(--muted);
        font-size: 0.76rem;
        font-weight: 750;
        text-transform: uppercase;
    }
    .contract-table td:first-child { color: var(--muted); font-weight: 650; }
    .contract-table td:last-child { color: var(--ink); }
    .contract-table tbody tr:hover { background: var(--surface-subtle); }
    .contract-table tbody tr:last-child td { border-bottom: 0; }
    .export-label { color: var(--muted); margin: 0.65rem 0 0.35rem; }
    section[data-testid="stMain"]::-webkit-scrollbar,
    [data-testid="stMainBlockContainer"]::-webkit-scrollbar {
        width: 10px;
    }
    section[data-testid="stMain"]::-webkit-scrollbar-thumb,
    [data-testid="stMainBlockContainer"]::-webkit-scrollbar-thumb {
        background: var(--scroll-thumb);
        border-radius: 5px;
    }
    section[data-testid="stMain"]::-webkit-scrollbar-track,
    [data-testid="stMainBlockContainer"]::-webkit-scrollbar-track {
        background: transparent;
    }

    @media (max-width: 1200px) {
        .st-key-global-reset-btn .stButton button {
            font-size: 0;
            left: auto;
            min-height: 2.35rem;
            padding: 0;
            right: calc(50% + min(22rem, calc((100vw - 26rem) / 2)) + 0.75rem);
            width: 2.35rem !important;
        }
        .st-key-global-reset-btn .stButton button::before {
            content: "\\21ba";
            display: block;
            font-family: "Segoe UI Symbol", sans-serif;
            font-size: 1rem;
            line-height: 1;
        }
    }

    @media (max-width: 900px) {
        section[data-testid="stSidebar"], section[data-testid="stSidebar"] > div {
            width: min(460px, 94vw) !important;
        }
        section[data-testid="stSidebar"] [data-testid="stHorizontalBlock"]:has(.st-key-scenario-action-rail) {
            gap: 0.5rem !important;
        }
        .block-container { padding: 4.5rem 5.4rem 3rem 1rem; }
        .model-item { border-right: 0; }
        /* 16px inputs stop iOS Safari auto-zooming into the field on focus
           (below 16px it zooms; this also covers phones in landscape). */
        section[data-testid="stSidebar"] input,
        .st-key-cost-panel input {
            font-size: 16px !important;
        }
        /* Main columns stack whenever the main area is narrow (Streamlit only
           wraps them on the narrowest screens, which leaves landscape phones
           and small tablets with squashed multi-column rows), and the metric
           row re-grids two-per-line. Streamlit's mobile column min-width
           (calc(100% - 24px)) would keep each metric card full-width, so it
           is zeroed here; the gap is pinned so the 50% basis arithmetic is
           deterministic. */
        [data-testid="stMainBlockContainer"] [data-testid="stHorizontalBlock"] {
            flex-wrap: wrap !important;
        }
        [data-testid="stMainBlockContainer"] [data-testid="stHorizontalBlock"]
        > [data-testid="stColumn"] {
            flex: 1 1 100% !important;
        }
        [data-testid="stMainBlockContainer"] [data-testid="stHorizontalBlock"]:has([data-testid="stMetric"]) {
            flex-wrap: wrap !important;
            gap: 0.75rem !important;
        }
        [data-testid="stMainBlockContainer"] [data-testid="stHorizontalBlock"]:has([data-testid="stMetric"])
        > [data-testid="stColumn"] {
            flex: 1 1 calc(50% - 0.375rem) !important;
            min-width: 0 !important;
            max-width: calc(50% - 0.375rem) !important;
        }
    }
    @media (max-width: 640px) {
        /* Phones: the sidebar becomes a full-width overlay sheet. Streamlit's
           responsive columns wrap the input column and the Calculate rail
           onto separate rows, which pushes the rail below the non-scrolling
           sidebar content and leaves Calculate unreachable. Force the rail's
           row into a column instead: the inputs scroll in the space above,
           and Calculate becomes a compact horizontal bar at the bottom — a
           fixed-height flex sibling, the same pattern as the desktop rail.
           flex: 0 0 auto defeats Streamlit's flex: 1 1 0% (flex-basis 0
           overrides the height), max-width: none beats its 90vw mobile cap,
           and 100dvh tracks the mobile URL bar better than 100vh. */
        section[data-testid="stSidebar"], section[data-testid="stSidebar"] > div {
            width: 100vw !important;
        }
        section[data-testid="stSidebar"] {
            max-width: none !important;
        }
        section[data-testid="stSidebar"] [data-testid="stHorizontalBlock"]:has(.st-key-scenario-action-rail) {
            flex: 0 0 auto !important;
            flex-direction: column !important;
            flex-wrap: nowrap !important;
            height: calc(100vh - 82px);
            height: calc(100dvh - 82px);
            gap: 0.6rem !important;
        }
        section[data-testid="stSidebar"] [data-testid="stHorizontalBlock"]:has(.st-key-scenario-action-rail)
        > [data-testid="stColumn"]:first-child {
            flex: 1 1 auto !important;
            height: auto !important;
            min-height: 0;
            overflow-y: auto;
            padding: 0 0 1rem;
        }
        section[data-testid="stSidebar"] .st-key-scenario-action-rail {
            flex: 0 0 auto !important;
            height: auto !important;
            min-height: 0 !important;
        }
        section[data-testid="stSidebar"] .st-key-scenario-action-rail
        [data-testid="stElementContainer"] {
            position: static;
            inset: auto;
        }
        section[data-testid="stSidebar"] .st-key-scenario-action-rail .stButton button {
            align-items: center;
            flex-direction: row;
            gap: 0.65rem;
            inset: auto;
            justify-content: center;
            padding: 0.8rem 1rem;
            position: static;
        }
        section[data-testid="stSidebar"] .st-key-scenario-action-rail .stButton button p {
            font-size: 1.1rem;
            white-space: nowrap;
            writing-mode: horizontal-tb;
            transform: none;
        }
        /* Streamlit collapses the sidebar with a fixed -300px translate; with
           this 100vw sheet that leaves a 90px sliver painted over the main
           content at z-index 999991. Hide the collapsed section outright (the
           expand control lives in the header, not the sidebar). */
        section[data-testid="stSidebar"][aria-expanded="false"] {
            visibility: hidden;
            transform: translateX(-100vw) !important;
        }
        /* Touch devices have no hover: keep the sidebar's close control
           visible so the sheet can always be dismissed. */
        section[data-testid="stSidebar"] [data-testid="stSidebarCollapseButton"] {
            visibility: visible;
        }
        /* The expand control is a bare 28px icon on the header — too easy to
           miss on a phone. Give it a labeled pill so the scenario builder is
           discoverable. */
        [data-testid="stExpandSidebarButton"] {
            align-items: center;
            background: var(--surface);
            border: 1px solid var(--line);
            border-radius: 999px;
            box-shadow: var(--shadow-sm);
            display: flex;
            gap: 0.45rem;
            min-height: 2.5rem;
            min-width: 2.5rem;
            padding: 0 0.85rem;
        }
        [data-testid="stExpandSidebarButton"]::after {
            content: "Edit inputs";
            color: var(--ink);
            font-size: 0.9rem;
            font-weight: 650;
        }

        /* Guards for real devices: no sideways scroll (an overflowing widget
           makes iOS/Android scroll the page right on focus, clipping the
           left edge), inputs that can shrink instead of overflowing the
           sheet, and metric values that wrap instead of clipping. */
        html, body { overflow-x: hidden; }
        section[data-testid="stSidebar"] [data-testid="stNumberInputContainer"],
        section[data-testid="stSidebar"] [data-testid="stSelectbox"],
        section[data-testid="stSidebar"] input {
            min-width: 0;
        }
        [data-testid="stMetricValue"], [data-testid="stMetricValue"] p {
            white-space: normal !important;
            overflow-wrap: anywhere;
        }
        /* Phones: the cost rail becomes a horizontal tab strip pinned to the
           bottom edge, and the panel opens upward as a bottom sheet. Extra
           bottom padding keeps the page content clear of the strip. */
        .st-key-cost-rail {
            bottom: 0.6rem;
            left: 0.6rem;
            right: 0.6rem;
            top: auto;
            width: auto;
        }
        .st-key-cost-rail [data-testid="stVerticalBlock"] {
            flex-direction: row !important;
            gap: 0.35rem !important;
        }
        .st-key-cost-rail [data-testid="stElementContainer"] {
            flex: 1;
            width: auto;
        }
        .st-key-cost-rail [role="radiogroup"] {
            flex-direction: row !important;
            flex-wrap: nowrap !important;
            gap: 0.3rem !important;
        }
        .st-key-cost-rail [role="radiogroup"] label {
            border-radius: 999px;
            border-right: 1px solid var(--line);
            box-shadow: var(--shadow-sm);
            flex: 1 1 0;
            font-size: 0.62rem;
            height: 2.6rem;
            min-height: 2.6rem;
            min-width: 0;
            padding: 0 0.3rem;
            transform: none;
            width: auto;
            writing-mode: horizontal-tb;
        }
        .st-key-cost-rail [role="radiogroup"] label:has(input[value="0"]) {
            flex: 0 0 2.4rem;
            font-size: 0.85rem;
            padding: 0 0.5rem;
        }
        .st-key-cost-rail [role="radiogroup"] label:has(input[value="0"])::after {
            content: none;
            margin-left: 0;
        }
        body:has(.st-key-cost_model_radio input[value="1"]:checked) .st-key-cost-rail label:has(input[value="1"]),
        body:has(.st-key-cost_model_radio input[value="2"]:checked) .st-key-cost-rail label:has(input[value="2"]),
        body:has(.st-key-cost_model_radio input[value="3"]:checked) .st-key-cost-rail label:has(input[value="3"]),
        body:has(.st-key-cost_model_radio input[value="4"]:checked) .st-key-cost-rail label:has(input[value="4"]),
        body:has(.st-key-cost_model_radio input[value="5"]:checked) .st-key-cost-rail label:has(input[value="5"]) {
            background: linear-gradient(168deg, var(--action-bg) 0%, var(--action-bg-hover) 130%);
            border-color: var(--action-bg);
            border-left: 1px solid var(--action-bg);
            color: var(--action-text);
            box-shadow: var(--shadow-sm);
            transform: translateY(-2px);
        }
        .st-key-cost-rail [role="radiogroup"] label::after {
            content: "";
            margin-left: auto;
            margin-top: 0;
        }
        .st-key-cost-rail [role="radiogroup"] label:active {
            transform: translateY(1px) scale(0.985);
        }
        .st-key-cost-panel {
            bottom: 4.6rem;
            left: 0.6rem;
            max-height: 55vh;
            right: 0.6rem;
            top: auto;
            transform: translateY(1.2rem) scale(0.99);
            width: auto;
        }
        body:has(.st-key-cost_model_radio input[value="1"]:checked) .st-key-cost-panel,
        body:has(.st-key-cost_model_radio input[value="2"]:checked) .st-key-cost-panel,
        body:has(.st-key-cost_model_radio input[value="3"]:checked) .st-key-cost-panel,
        body:has(.st-key-cost_model_radio input[value="4"]:checked) .st-key-cost-panel,
        body:has(.st-key-cost_model_radio input[value="5"]:checked) .st-key-cost-panel {
            transform: none;
        }
        /* Phones keep the sheet's bottom Calculate bar and hide the header
           button (the header is too crowded at phone width). The desktop rail
           layout (an absolute button filling a full-height rail) is reverted
           to an in-flow compact bar so nothing overflows the viewport. */
        section[data-testid="stSidebar"] .st-key-scenario-action-rail {
            display: block;
            flex: 0 0 auto !important;
            height: auto !important;
            min-height: 0 !important;
        }
        section[data-testid="stSidebar"] .st-key-scenario-action-rail
        [data-testid="stElementContainer"] {
            position: static;
            inset: auto;
        }
        section[data-testid="stSidebar"] .st-key-scenario-action-rail [data-testid="stFormSubmitButton"] button {
            align-items: center;
            flex-direction: row;
            gap: 0.65rem;
            inset: auto;
            justify-content: center;
            padding: 0.8rem 1rem;
            position: static;
            width: 100% !important;
        }
        section[data-testid="stSidebar"] .st-key-scenario-action-rail
        [data-testid="stFormSubmitButton"] button::before,
        section[data-testid="stSidebar"] .st-key-scenario-action-rail
        [data-testid="stFormSubmitButton"] button::after {
            content: "\\2192";
            display: block;
            flex: 0 0 auto;
            font-family: "Segoe UI Symbol", sans-serif;
            font-size: 1.2rem;
            line-height: 1;
        }
        section[data-testid="stSidebar"] .st-key-scenario-action-rail [data-testid="stFormSubmitButton"] button p {
            font-size: 1.1rem;
            white-space: nowrap;
            writing-mode: horizontal-tb;
            transform: none;
        }
        .st-key-calculate-anchor {
            display: none;
        }
        .st-key-global-reset-btn {
            display: none;
        }
        .st-key-calculate_panel {
            background: var(--surface);
            bottom: -1.2rem;
            display: block;
            padding-top: 0.45rem;
            position: sticky;
            z-index: 2;
        }
        .st-key-calculate_panel [data-testid="stFormSubmitButton"] button {
            background: var(--action-bg);
            border-color: var(--action-bg);
            border-radius: 10px;
            color: var(--action-text);
            font-weight: 700;
            min-height: 2.75rem;
            width: 100% !important;
        }
        .block-container { padding: 4.25rem 0.75rem 6.5rem; }
        h1 { font-size: 1.8rem; }
        .page-deck { font-size: 0.92rem; }
        [data-testid="stMetric"] { min-height: 96px; padding: 0.8rem 0.75rem; }
        [data-testid="stMetricValue"] { font-size: 1.18rem; }
        .contract-card { overflow-x: auto; }
        div[data-testid="stTabs"] div[role="tablist"] { gap: 1.1rem; }
        .scenario-bar { align-items: flex-start; flex-direction: column; gap: 0.25rem; }
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# Theme: an explicit choice (URL parameter or the sidebar toggle) wins; on a
# first visit Python starts from light and the component below redirects to
# dark when the system prefers it. The rendered marker drives the CSS via
# html:has(...), so every rerun reflects the resolved theme immediately.
theme = str(st.query_params.get("theme", ""))
if theme not in ("light", "dark"):
    theme = st.session_state.get("theme", "light")
    if theme not in ("light", "dark"):
        theme = "light"
st.session_state["theme"] = theme

st.markdown(
    f'<div class="theme-marker theme-{theme}" hidden></div>',
    unsafe_allow_html=True,
)


def _apply_system_dark() -> None:
    st.session_state["theme"] = "dark"
    st.query_params["theme"] = "dark"
    st.rerun()


# Hidden button the theme-sync iframe clicks programmatically when the system
# prefers dark and no explicit theme is set, so the server learns the
# preference (its iframe is sandboxed against navigating the top window).
st.button(
    label=None,
    key="theme_auto_dark",
    on_click=_apply_system_dark,
    icon=":material/dark_mode:",
)

_THEME_SYNC_JS = """
<script>
(() => {
  const params = new URLSearchParams(parent.location.search);
  if (params.get("theme")) return;
  if (!parent.matchMedia("(prefers-color-scheme: dark)").matches) return;
  // Apply the dark classes to the theme marker immediately: the whole CSS
  // dark palette hangs off html:has(.theme-marker.theme-dark), so flipping
  // the class re-themes the page in the same frame instead of waiting for
  // the server rerun (which is what made the tabs flash light for seconds).
  const tryApply = () => {
    const marker = parent.document.querySelector(".theme-marker");
    if (!marker) return false;
    marker.classList.remove("theme-light");
    marker.classList.add("theme-dark");
    return true;
  };
  const tryClick = () => {
    // Still click the hidden button so the SERVER also learns the theme and
    // the rerun re-renders charts and markers in dark consistently.
    const btn = parent.document.querySelector(".st-key-theme_auto_dark button");
    if (!btn) return false;
    btn.click();
    return true;
  };
  tryApply();
  if (tryClick()) return;
  let tries = 0;
  const timer = setInterval(() => {
    tries += 1;
    if (tryClick() || tries > 100) clearInterval(timer);
  }, 150);
})();
</script>
"""

st.iframe(_THEME_SYNC_JS, width=1, height=1)

# Click-away close for the cost panel: a one-shot iframe attaches a document
# listener on the parent page that selects the rail's ✕ radio option whenever
# a click lands outside the rail and the panel. The flag guards against
# duplicate listeners (the iframe is re-created on every rerun).
_COST_CLICK_AWAY_JS = """
<script>
(() => {
  if (parent.document.__costClickAway) return;
  parent.document.__costClickAway = true;
  let downInside = false;
  const inside = (target) => {
    const rail = parent.document.querySelector(".st-key-cost-rail");
    const panel = parent.document.querySelector(".st-key-cost-panel");
    return (rail && rail.contains(target)) || (panel && panel.contains(target));
  };
  // A drag that STARTS inside the rail or panel must never close the panel,
  // even if it is released outside (the click then fires on the common
  // ancestor, which would look like an outside click).
  parent.document.addEventListener("pointerdown", (event) => {
    downInside = inside(event.target);
  }, true);
  parent.document.addEventListener("click", (event) => {
    const closePanel = () => {
      const closeLabel = parent.document.querySelector(
        '.st-key-cost-rail label:has(input[value="0"])'
      );
      // Deferred: a click synthesized inside another event's capture phase
      // is swallowed, so close after the original click has fully completed.
      if (closeLabel) setTimeout(() => closeLabel.click(), 0);
    };
    const targetInside = inside(event.target);
    if (downInside) {
      // Released outside after pressing inside: keep the panel open.
      if (!targetInside) return;
      // Pressed and released on the already-selected model tab: toggle close.
      const label = event.target.closest(".st-key-cost-rail [role='radiogroup'] label");
      if (label) {
        const input = label.querySelector("input");
        if (input && input.checked && input.value !== "0") closePanel();
      }
      return;
    }
    if (targetInside) return;
    closePanel();
  }, true);
})();
</script>
"""
st.iframe(_COST_CLICK_AWAY_JS, width=1, height=1)

# On phones the sidebar is a full-screen sheet, so after Calculate the fresh
# results render hidden behind it. A new iframe (rendered only on the post-
# submit rerun) marks the sidebar collapsed in Streamlit's localStorage and
# reloads the page: the reload is a same-document navigation, which the
# iframe sandbox permits (redirects elsewhere are blocked), and the reloaded
# app re-runs from the URL params the submit just wrote, so the committed
# results are what render. Desktop skips this: its sidebar sits beside the
# content, so the results are already visible.
_SIDEBAR_CLOSE_JS = """
<script>
(() => {
  if (parent.innerWidth > 640) return;
  const key = Object.keys(parent.localStorage).find((k) => k.startsWith("stSidebarCollapsed"))
    || "stSidebarCollapsed-";
  parent.localStorage.setItem(key, "true");
  parent.location.reload();
})();
</script>
"""

# Theme toggle, pinned to the page header next to Streamlit's toolbar menu.
with st.container(key="theme-toggle-anchor"):
    if st.button(
        "Dark" if theme == "light" else "Light",
        icon=":material/dark_mode:" if theme == "light" else ":material/light_mode:",
        key="theme_toggle",
    ):
        new_theme = "dark" if theme == "light" else "light"
        st.session_state["theme"] = new_theme
        st.query_params["theme"] = new_theme
        st.rerun()

# Cost-assumption models for the right-edge rail. Each key is a widget prefix
# (matching the scenario parameter names), and the label is the tab handle.
COST_MODELS = [
    ("dna", "DNA"),
    ("amazon", "Amazon"),
    ("azure", "Azure"),
    ("tape", "Tape"),
    ("custom", "Custom"),
]

# Widget keys per cost model, used by the per-model "Reset to defaults"
# buttons: only the listed keys are restored to the paper baseline.
MODEL_WIDGET_KEYS = {
    "dna": [
        "dna_cost_base_year", "dna_synthesis_cost", "dna_sequencing_cost",
        "synthesis_decline", "sequencing_decline", "dna_durability",
    ],
    "amazon": [
        "amazon_base_year", "amazon_put_per_1000", "amazon_restore_per_1000",
        "amazon_retrieval_per_tb", "amazon_storage_per_tb_month", "amazon_decline",
    ],
    "azure": [
        "azure_base_year", "azure_write_per_1000", "azure_read_per_1000",
        "azure_retrieval_per_tb", "azure_storage_per_tb_month", "azure_decline",
    ],
    "tape": [
        "tape_base_year", "tape_durability", "tape_media_per_tb",
        "tape_hardware_per_tb", "tape_energy_per_tb_year", "tape_media_decline",
        "tape_hardware_decline", "tape_energy_decline",
    ],
    "custom": [
        "custom_name", "custom_base_year", "custom_write_tb", "custom_write_asset",
        "custom_storage_tb_year", "custom_retrieval_tb", "custom_retrieval_asset",
        "custom_decline", "custom_replacement",
    ],
}



SIDEBAR_WIDGET_KEYS = [
    "archive_value",
    "archive_unit",
    "asset_value",
    "asset_unit",
    "retrieval",
    "start_year_widget",
    "horizon",
    "tech_dna",
    "tech_amazon",
    "tech_azure",
    "tech_tape",
    "tech_custom",
    "discount",
    "projection_end",
    "log_scale",
]


WIDGET_KEYS = [
    "archive_value",
    "archive_unit",
    "asset_value",
    "asset_unit",
    "retrieval",
    "start_year_widget",
    "horizon",
    "tech_dna",
    "tech_amazon",
    "tech_azure",
    "tech_tape",
    "tech_custom",
    "discount",
    "dna_cost_base_year",
    "dna_synthesis_cost",
    "dna_sequencing_cost",
    "synthesis_decline",
    "sequencing_decline",
    "amazon_base_year",
    "amazon_put_per_1000",
    "amazon_restore_per_1000",
    "amazon_retrieval_per_tb",
    "amazon_storage_per_tb_month",
    "amazon_decline",
    "azure_base_year",
    "azure_write_per_1000",
    "azure_read_per_1000",
    "azure_retrieval_per_tb",
    "azure_storage_per_tb_month",
    "azure_decline",
    "tape_base_year",
    "tape_media_per_tb",
    "tape_hardware_per_tb",
    "tape_energy_per_tb_year",
    "tape_media_decline",
    "tape_hardware_decline",
    "tape_energy_decline",
    "dna_durability",
    "tape_durability",
    "custom_name",
    "custom_base_year",
    "custom_write_tb",
    "custom_write_asset",
    "custom_storage_tb_year",
    "custom_retrieval_tb",
    "custom_retrieval_asset",
    "custom_decline",
    "custom_replacement",
    "projection_end",
    "log_scale",
]


def _default_projection_end(scenario: Scenario, query: dict[str, str]) -> int:
    return min(2500, max(scenario.start_year + 25, int(query.get("projection_end", 2350))))


def _query_mapping() -> dict[str, str]:
    return {key: str(value) for key, value in st.query_params.items()}


def _initial_scenario() -> Scenario:
    try:
        return Scenario.from_mapping(_query_mapping())
    except (TypeError, ValueError):
        st.warning("The shared URL contained invalid values. The paper baseline has been loaded instead.")
        return Scenario()


def _display_units(scenario: Scenario) -> tuple[float, str, float, str]:
    if scenario.archive_size_tb >= 1_000_000:
        archive_value, archive_unit = scenario.archive_size_tb / 1_000_000, "EB"
    elif scenario.archive_size_tb >= 1_000:
        archive_value, archive_unit = scenario.archive_size_tb / 1_000, "PB"
    else:
        archive_value, archive_unit = scenario.archive_size_tb, "TB"
    if scenario.average_asset_size_mb >= 1000:
        asset_value, asset_unit = scenario.average_asset_size_mb / 1000, "GB"
    else:
        asset_value, asset_unit = scenario.average_asset_size_mb, "MB"
    return archive_value, archive_unit, asset_value, asset_unit


def _widget_state_from_scenario(
    scenario: Scenario, query: dict[str, str] | None = None
) -> dict[str, bool | float | int | str]:
    archive_value, archive_unit, asset_value, asset_unit = _display_units(scenario)
    query = query or {}
    return {
        "archive_value": archive_value,
        "archive_unit": archive_unit,
        "asset_value": asset_value,
        "asset_unit": asset_unit,
        "retrieval": scenario.annual_retrieval_percent,
        "start_year_widget": scenario.start_year,
        "horizon": scenario.horizon_years,
        "tech_dna": "DNA" in scenario.technologies,
        "tech_amazon": "Amazon Deep Archive" in scenario.technologies,
        "tech_azure": "Azure Blob Archive" in scenario.technologies,
        "tech_tape": "Tape On-premise" in scenario.technologies,
        "tech_custom": "Custom storage" in scenario.technologies,
        "discount": scenario.discount_rate_percent,
        "dna_cost_base_year": scenario.dna_cost_base_year,
        "dna_synthesis_cost": scenario.dna_synthesis_cost_per_mb,
        "dna_sequencing_cost": scenario.dna_sequencing_cost_per_mb,
        "synthesis_decline": scenario.synthesis_decline_percent,
        "sequencing_decline": scenario.sequencing_decline_percent,
        "amazon_base_year": scenario.amazon_price_base_year,
        "amazon_put_per_1000": scenario.amazon_put_usd_per_request * 1_000,
        "amazon_restore_per_1000": scenario.amazon_bulk_restore_usd_per_request * 1_000,
        "amazon_retrieval_per_tb": scenario.amazon_bulk_retrieval_usd_per_mb * 1_000_000,
        "amazon_storage_per_tb_month": scenario.amazon_storage_usd_per_mb_month * 1_000_000,
        "amazon_decline": scenario.amazon_decline_percent,
        "azure_base_year": scenario.azure_price_base_year,
        "azure_write_per_1000": scenario.azure_write_usd_per_request * 1_000,
        "azure_read_per_1000": scenario.azure_read_usd_per_request * 1_000,
        "azure_retrieval_per_tb": scenario.azure_retrieval_usd_per_mb * 1_000_000,
        "azure_storage_per_tb_month": scenario.azure_storage_usd_per_mb_month * 1_000_000,
        "azure_decline": scenario.azure_decline_percent,
        "tape_base_year": scenario.tape_price_base_year,
        "tape_media_per_tb": scenario.tape_media_usd_per_tb,
        "tape_hardware_per_tb": scenario.tape_hardware_usd_per_tb,
        "tape_energy_per_tb_year": scenario.tape_energy_usd_per_tb_year,
        "tape_media_decline": scenario.tape_media_decline_percent,
        "tape_hardware_decline": scenario.tape_hardware_decline_percent,
        "tape_energy_decline": scenario.tape_energy_decline_percent,
        "dna_durability": scenario.dna_durability_years,
        "tape_durability": scenario.tape_durability_years,
        "custom_name": scenario.custom_storage_name,
        "custom_base_year": scenario.custom_cost_base_year,
        "custom_write_tb": scenario.custom_write_cost_per_tb,
        "custom_write_asset": scenario.custom_write_cost_per_asset,
        "custom_storage_tb_year": scenario.custom_storage_cost_per_tb_year,
        "custom_retrieval_tb": scenario.custom_retrieval_cost_per_tb,
        "custom_retrieval_asset": scenario.custom_retrieval_cost_per_asset,
        "custom_decline": scenario.custom_decline_percent,
        "custom_replacement": scenario.custom_replacement_years,
        "projection_end": _default_projection_end(scenario, query),
        "log_scale": query.get("log_scale", "True").lower() == "true",
    }


def _snapshot_widgets() -> dict[str, bool | float | int | str]:
    return {key: st.session_state[key] for key in WIDGET_KEYS}


def _scenario_from_widgets(widgets: dict[str, bool | float | int | str]) -> Scenario:
    archive_multiplier = {"TB": 1, "PB": 1000, "EB": 1_000_000}[widgets["archive_unit"]]
    asset_multiplier = {"MB": 1, "GB": 1000}[widgets["asset_unit"]]
    technologies = tuple(
        technology
        for selected, technology in (
            (widgets["tech_dna"], "DNA"),
            (widgets["tech_amazon"], "Amazon Deep Archive"),
            (widgets["tech_azure"], "Azure Blob Archive"),
            (widgets["tech_tape"], "Tape On-premise"),
            (widgets["tech_custom"], "Custom storage"),
        )
        if selected
    )
    return Scenario(
        archive_size_tb=widgets["archive_value"] * archive_multiplier,
        average_asset_size_mb=widgets["asset_value"] * asset_multiplier,
        annual_retrieval_percent=widgets["retrieval"],
        start_year=int(widgets["start_year_widget"]),
        horizon_years=int(widgets["horizon"]),
        discount_rate_percent=widgets["discount"],
        dna_cost_base_year=int(widgets["dna_cost_base_year"]),
        dna_synthesis_cost_per_mb=widgets["dna_synthesis_cost"],
        dna_sequencing_cost_per_mb=widgets["dna_sequencing_cost"],
        synthesis_decline_percent=widgets["synthesis_decline"],
        sequencing_decline_percent=widgets["sequencing_decline"],
        amazon_price_base_year=int(widgets["amazon_base_year"]),
        amazon_put_usd_per_request=widgets["amazon_put_per_1000"] / 1_000,
        amazon_bulk_restore_usd_per_request=widgets["amazon_restore_per_1000"] / 1_000,
        amazon_bulk_retrieval_usd_per_mb=widgets["amazon_retrieval_per_tb"] / 1_000_000,
        amazon_storage_usd_per_mb_month=widgets["amazon_storage_per_tb_month"] / 1_000_000,
        amazon_decline_percent=widgets["amazon_decline"],
        azure_price_base_year=int(widgets["azure_base_year"]),
        azure_write_usd_per_request=widgets["azure_write_per_1000"] / 1_000,
        azure_read_usd_per_request=widgets["azure_read_per_1000"] / 1_000,
        azure_retrieval_usd_per_mb=widgets["azure_retrieval_per_tb"] / 1_000_000,
        azure_storage_usd_per_mb_month=widgets["azure_storage_per_tb_month"] / 1_000_000,
        azure_decline_percent=widgets["azure_decline"],
        tape_price_base_year=int(widgets["tape_base_year"]),
        tape_media_usd_per_tb=widgets["tape_media_per_tb"],
        tape_hardware_usd_per_tb=widgets["tape_hardware_per_tb"],
        tape_energy_usd_per_tb_year=widgets["tape_energy_per_tb_year"],
        tape_media_decline_percent=widgets["tape_media_decline"],
        tape_hardware_decline_percent=widgets["tape_hardware_decline"],
        tape_energy_decline_percent=widgets["tape_energy_decline"],
        dna_durability_years=int(widgets["dna_durability"]),
        tape_durability_years=int(widgets["tape_durability"]),
        custom_storage_name=widgets["custom_name"],
        custom_cost_base_year=int(widgets["custom_base_year"]),
        custom_write_cost_per_tb=widgets["custom_write_tb"],
        custom_write_cost_per_asset=widgets["custom_write_asset"],
        custom_storage_cost_per_tb_year=widgets["custom_storage_tb_year"],
        custom_retrieval_cost_per_tb=widgets["custom_retrieval_tb"],
        custom_retrieval_cost_per_asset=widgets["custom_retrieval_asset"],
        custom_decline_percent=widgets["custom_decline"],
        custom_replacement_years=int(widgets["custom_replacement"]),
        technologies=technologies,
    )


def _money(value: float) -> str:
    absolute = abs(value)
    for threshold, suffix in ((1e12, "T"), (1e9, "B"), (1e6, "M"), (1e3, "K")):
        if absolute >= threshold:
            return f"${value / threshold:,.2f}{suffix}"
    if absolute >= 1:
        return f"${value:,.2f}"
    return f"${value:,.3g}"


def _quantity(value: float) -> str:
    absolute = abs(value)
    for threshold, suffix in ((1e15, "Q"), (1e12, "T"), (1e9, "B"), (1e6, "M"), (1e3, "K")):
        if absolute >= threshold:
            return f"{value / threshold:,.2f}{suffix}"
    return f"{value:,.0f}"


def _section_intro(title: str, description: str) -> None:
    st.markdown(
        f'<div class="tab-intro"><h2>{title}</h2><p>{description}</p></div>',
        unsafe_allow_html=True,
    )


def _themed_table(columns: list[str], rows: list[list[str]]) -> None:
    header = "".join(f"<th>{html.escape(column)}</th>" for column in columns)
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(cell))}</td>" for cell in row) + "</tr>"
        for row in rows
    )
    st.markdown(
        f'<div class="contract-card"><table class="contract-table">'
        f"<thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div>",
        unsafe_allow_html=True,
    )


def _chart_downloads(
    key: str,
    csv_data: bytes,
    filename_base: str,
    chart_key: str,
) -> None:
    st.markdown('<div class="export-label">Export chart</div>', unsafe_allow_html=True)
    columns = st.columns(3)
    columns[0].download_button(
        "CSV",
        csv_data,
        f"{filename_base}.csv",
        "text/csv",
        key=f"download_{key}_csv",
        icon=":material/download:",
        width="stretch",
        help="Download the data shown in this graph.",
    )
    columns[1].button(
        "PNG",
        key=f"download_{key}_png",
        icon=":material/download:",
        width="stretch",
        help="Download this graph as a high-resolution PNG image.",
    )
    columns[2].button(
        "SVG",
        key=f"download_{key}_svg",
        icon=":material/download:",
        width="stretch",
        help="Download this graph as an editable vector image.",
    )
    _bind_live_plotly_image_downloads(key, chart_key, filename_base)


def _bind_live_plotly_image_downloads(key: str, chart_key: str, filename_base: str) -> None:
    config_json = json.dumps(
        {
            "chartKey": chart_key,
            "filename": filename_base,
            "buttons": {
                "png": f"download_{key}_png",
                "svg": f"download_{key}_svg",
            },
        }
    )
    st.iframe(
        f"""
        <script>
        (() => {{
          const config = {config_json};
          const findByKey = (key) => parent.document.querySelector(`.st-key-${{key}}`);
          const findGraph = () => {{
            const root = findByKey(config.chartKey);
            return root ? root.querySelector(".js-plotly-plot") : null;
          }};
          const plotlyApi = () => parent.Plotly || parent.document.defaultView.Plotly;
          const bind = () => {{
            const graph = findGraph();
            const Plotly = plotlyApi();
            let bound = 0;
            for (const [format, buttonKey] of Object.entries(config.buttons)) {{
              const buttonRoot = findByKey(buttonKey);
              const button = buttonRoot ? buttonRoot.querySelector("button") : null;
              if (!button || button.dataset.livePlotlyExport === config.chartKey + format) {{
                if (button) bound += 1;
                continue;
              }}
              button.dataset.livePlotlyExport = config.chartKey + format;
              button.addEventListener("click", (event) => {{
                const liveGraph = findGraph();
                const livePlotly = plotlyApi();
                if (!liveGraph || !livePlotly || !livePlotly.downloadImage) return;
                event.preventDefault();
                event.stopPropagation();
                event.stopImmediatePropagation();
                const bounds = liveGraph.getBoundingClientRect();
                const width = Math.max(320, Math.round(bounds.width || liveGraph._fullLayout?.width || 1200));
                const height = Math.max(240, Math.round(bounds.height || liveGraph._fullLayout?.height || 500));
                livePlotly.downloadImage(liveGraph, {{
                  format,
                  filename: config.filename,
                  width,
                  height,
                  scale: format === "png" ? 2 : 1,
                }});
              }}, true);
              bound += 1;
            }}
            return Boolean(graph && Plotly && Plotly.downloadImage && bound === 2);
          }};
          if (bind()) return;
          let tries = 0;
          const timer = setInterval(() => {{
            tries += 1;
            if (bind() || tries > 100) clearInterval(timer);
          }}, 100);
        }})();
        </script>
        """,
        width=1,
        height=1,
    )


def _plot_config(filename: str) -> dict:
    return {
        "displaylogo": False,
        "toImageButtonOptions": {"format": "png", "filename": filename, "scale": 2},
    }


def _plotly_chart_with_placeholder(figure, *, key: str, filename: str) -> None:
    slot = st.empty()
    height = int(figure.layout.height or 500)
    with slot.container():
        st.markdown(
            f'<div class="chart-loading" style="height: {height}px;"></div>',
            unsafe_allow_html=True,
        )
    slot.plotly_chart(
        figure,
        key=key,
        width="stretch",
        config=_plot_config(filename),
    )


@st.cache_data(show_spinner=False, max_entries=32)
def _cached_simulation(scenario: Scenario):
    return simulate_scenario(scenario)


@st.cache_data(show_spinner=False, max_entries=32)
def _cached_projection(scenario: Scenario, final_year: int) -> pd.DataFrame:
    return simulate_start_years(scenario, final_year)


@st.cache_data(show_spinner=False, max_entries=32)
def _cached_dna_costs(scenario: Scenario, final_year: int) -> pd.DataFrame:
    return simulate_dna_unit_costs(scenario, final_year)


initial = _initial_scenario()
query = _query_mapping()
initial_widgets = _widget_state_from_scenario(initial, query)
baseline_widgets = _widget_state_from_scenario(Scenario())


def _reset_widget_keys(widget_keys: list[str]) -> None:
    for widget_key in widget_keys:
        st.session_state[widget_key] = baseline_widgets[widget_key]


for widget_key, default_value in initial_widgets.items():
    st.session_state.setdefault(widget_key, default_value)
# The graphs only follow the last calculated inputs: the initial load
# counts as the first calculation. Every widget lives inside one form, so
# edits never trigger reruns — the three Calculate buttons are form submits
# and are the only way to commit changes.
st.session_state.setdefault("committed_widgets", dict(initial_widgets))

# All inputs live in one form: Streamlit batches form widgets client-side and
# reruns the script only on a submit, so editing parameters never re-renders
# the buttons or the charts. The form block lives in the sidebar; the panel
# and header Calculate join it via the main dg's form data below.
with st.sidebar:
    with st.container(key="sidebar-reset-btn"):
        st.button(
            "Reset sidebar values to paper baseline",
            key="sidebar_reset",
            on_click=_reset_widget_keys,
            args=(SIDEBAR_WIDGET_KEYS,),
            width="stretch",
        )
    with st.form("scenario_form", border=False, enter_to_submit=False):
        input_column, action_column = st.columns([6, 1], gap="small")
        with input_column:
            st.markdown(
                """
                <div class="sidebar-kicker">Scenario builder</div>
                <div class="sidebar-title">Model inputs</div>
                <p class="sidebar-copy">Archive workload, time horizon, and technology selection. Cost assumptions live in the tabs on the right edge of the page.</p>
                """,
                unsafe_allow_html=True,
            )
            st.markdown('<div class="sidebar-section">Workload and time</div>', unsafe_allow_html=True)
            col_a, col_b = st.columns([2, 1])
            with col_a:
                archive_input = st.number_input(
                    "Archive size",
                    min_value=0.001,
                    key="archive_value",
                    help="Total logical data stored in the collection, before DNA coding, redundancy, or provider replication.",
                )
            with col_b:
                archive_unit_input = st.selectbox(
                    "Unit", ["TB", "PB", "EB"], key="archive_unit",
                    help="Decimal capacity unit: 1 PB = 1,000 TB and 1 EB = 1,000,000 TB.",
                )
    
            col_a, col_b = st.columns([2, 1])
            with col_a:
                asset_input = st.number_input(
                    "Average asset size",
                    min_value=0.001,
                    key="asset_value",
                    help="Average data object size. It determines object count.",
                )
            with col_b:
                asset_unit_input = st.selectbox(
                    "Unit ", ["MB", "GB"], key="asset_unit",
                    help="Unit used for the average size of one asset.",
                )
    
            time_col_a, time_col_b = st.columns(2)
            with time_col_a:
                start_year = st.number_input(
                    "Start year", min_value=2025, max_value=2500,
                    key="start_year_widget", help="Calendar year in which the archive is first written.",
                )
            with time_col_b:
                horizon = st.number_input(
                    "Retention (years)", min_value=1, max_value=10_000,
                    key="horizon",
                    help="Number of charged storage years, including the start year.",
                )
    
            finance_col_a, finance_col_b = st.columns(2)
            with finance_col_a:
                retrieval = st.number_input(
                    "Annual retrieval (%)", min_value=0.0, max_value=10_000.0,
                    step=0.25, key="retrieval",
                    help="Expected share of the logical archive retrieved each year. 1% means reading 10 TB per year from a 1 PB archive.",
                )
            with finance_col_b:
                discount = st.number_input(
                    "Discount rate (%)", min_value=0.0, max_value=99.0,
                    step=0.25, key="discount",
                    help="Real rate used to discount future payments to the storage start year. Use 0 for undiscounted costs.",
                )
    
            st.markdown('<div class="sidebar-section">Display</div>', unsafe_allow_html=True)
            chart_col_a, chart_col_b = st.columns([1.4, 1])
            with chart_col_a:
                projection_end = st.number_input(
                    "Outlook end year", min_value=int(start_year), max_value=2500,
                    key="projection_end",
                    help="Final archive start year included in the start-year outlook chart.",
                )
            with chart_col_b:
                st.markdown('<div class="log-scale-label">Log</div>', unsafe_allow_html=True)
                log_scale = st.toggle(
                    label="Log scale",
                    label_visibility="collapsed",
                    key="log_scale",
                    help="Recommended when technologies differ by several orders of magnitude.",
                )
    
            st.markdown('<div class="sidebar-section">Technologies</div>', unsafe_allow_html=True)
            with st.container(border=True):
                tech_col_a, tech_col_b = st.columns(2)
                with tech_col_a:
                    tech_dna = st.checkbox(
                        "DNA", key="tech_dna", help="Archival storage using DNA synthesis and sequencing.",
                    )
                    tech_amazon = st.checkbox(
                        "Amazon S3", key="tech_amazon", help="Amazon S3 Glacier Deep Archive.",
                    )
                    tech_tape = st.checkbox(
                        "Tape", key="tech_tape", help="On-premise tape storage with periodic media replacement.",
                    )
                with tech_col_b:
                    tech_azure = st.checkbox(
                        "Azure Blob", key="tech_azure", help="Microsoft Azure Blob Storage Archive tier.",
                    )
                    tech_custom = st.checkbox(
                        "Custom", key="tech_custom", help="User-defined storage cost model.",
                    )
    
            st.markdown('<div class="sidebar-section">Cost assumptions</div>', unsafe_allow_html=True)
            st.markdown(
                '<p class="sidebar-copy">Editable unit costs and decline rates for each storage model live in the tabs on the right edge of the page.</p>',
                unsafe_allow_html=True,
            )
        with action_column:
            with st.container(key="scenario-action-rail"):
                # The mobile bottom bar; hidden on desktop, where the wide
                # centered header button is the primary Calculate.
                calculate_scenario = st.form_submit_button(
                    "Calculate",
                    type="primary",
                    key="calculate_scenario",
                    width="stretch",
                )
    
# The panel and header Calculate join the sidebar's form. st.form cannot
# wrap both containers (a form is a single block), so the form id is attached
# to the main dg directly — the same mechanism the form block uses on itself.
# Widgets added through `st.foo` calls read this dg's form data, so the panel
# inputs batch with the sidebar's and only a Calculate submit reruns anything.
# Cost-assumption rail and panel: slim vertical model tabs on the right
# edge of the page (a bottom strip on phones). The tabs are pure-HTML radio
# labels — opening, closing, and switching never trigger a script rerun, so
# the buttons and charts stay untouched. All five model blocks stay mounted
# and CSS shows only the checked one, so edits survive closing and switching.
# Wide centered Calculate button in the top header: it stays visible while
# the sidebar or a cost panel is open, so it works after either kind of edit.
# (Hidden on phones, where the sheet bar and the panel bar cover the flows.)
with st.container(key="global-reset-btn"):
    st.button(
        "Reset all inputs to paper baseline",
        key="global_reset",
        on_click=_reset_widget_keys,
        args=(WIDGET_KEYS,),
    )

_main_dg._form_data = _FormData("scenario_form")


def _render_model_reset(model_key: str) -> None:
    with st.container(key=f"model-reset-{model_key}"):
        st.form_submit_button(
            "Reset to defaults",
            key=f"reset_{model_key}",
            on_click=_reset_widget_keys,
            args=(MODEL_WIDGET_KEYS[model_key],),
            icon=":material/refresh:",
        )


with st.container(key="calculate-anchor"):
    calculate_header = st.form_submit_button(
        "Calculate",
        key="calculate_header",
        type="primary",
        width="stretch",
    )

# The tabs are a Streamlit radio: the frontend manages its checked state
# instantly (no rerun \u2014 it is a form widget), the radio group enforces
# exclusivity natively, and the CSS reads the checked input's value to open
# the matching panel. "\u2715" is the closed state. (Raw HTML radios/details do
# not work: Streamlit's page scripts suppress native form-control activation.)
st.session_state.setdefault("cost_model_radio", "\u2715")

with st.container(key="cost-rail"):
    st.radio(
        "Cost model",
        options=["\u2715", "DNA", "Amazon", "Azure", "Tape", "Custom"],
        key="cost_model_radio",
        horizontal=False,
        label_visibility="collapsed",
    )

with st.container(key="cost-panel"):
    st.markdown(
        '<div class="cost-panel-header">'
        '<div class="sidebar-kicker">Cost assumptions</div>'
        '<div class="cost-panel-title cost-title-none">Cost assumptions</div>'
        '<div class="cost-panel-title cost-title-dna">DNA cost assumptions</div>'
        '<div class="cost-panel-title cost-title-amazon">Amazon Deep Archive assumptions</div>'
        '<div class="cost-panel-title cost-title-azure">Azure Blob Archive assumptions</div>'
        '<div class="cost-panel-title cost-title-tape">Tape on-premise assumptions</div>'
        '<div class="cost-panel-title cost-title-custom">Custom storage assumptions</div>'
        '<p class="cost-panel-deck">Editable unit costs, decline rates, and replacement cycles for each storage model. Changes apply when you press Calculate.</p>'
        '<div class="cost-panel-summary cost-summary-dna"><span>Synthesis</span><span>Sequencing</span><span>Durability</span></div>'
        '<div class="cost-panel-summary cost-summary-amazon"><span>Requests</span><span>Retrieval</span><span>Storage</span></div>'
        '<div class="cost-panel-summary cost-summary-azure"><span>Operations</span><span>Retrieval</span><span>Archive tier</span></div>'
        '<div class="cost-panel-summary cost-summary-tape"><span>Media</span><span>Hardware</span><span>Energy</span></div>'
        '<div class="cost-panel-summary cost-summary-custom"><span>Write</span><span>Store</span><span>Retrieve</span></div>'
        '</div>',
        unsafe_allow_html=True,
    )

    with st.container(key="cost_model_dna"):
        _render_model_reset("dna")
        dna_cost_base_year = st.number_input(
            "DNA cost base year", min_value=2000, max_value=2500,
            key="dna_cost_base_year",
            help="Year to which the editable synthesis and sequencing unit costs apply.",
        )
        dna_synthesis_cost = st.number_input(
            "Synthesis cost (USD/MB)", min_value=0.0,
            format="%.6f", key="dna_synthesis_cost",
            help="Cost in the DNA cost base year to synthesize enough bases for 1 MB of logical data, before redundancy and indexing overhead.",
        )
        dna_sequencing_cost = st.number_input(
            "Sequencing cost (USD/MB)", min_value=0.0,
            format="%.8f", key="dna_sequencing_cost",
            help="Cost in the DNA cost base year to sequence 1 MB of retrieved logical data.",
        )
        synthesis_decline = st.number_input(
            "Synthesis annual decline (%)", min_value=0.0, max_value=99.99,
            key="synthesis_decline",
            help="Percentage by which synthesis cost is assumed to fall each calendar year.",
        )
        sequencing_decline = st.number_input(
            "Sequencing annual decline (%)", min_value=0.0, max_value=99.99,
            key="sequencing_decline",
            help="Percentage by which sequencing cost is assumed to fall each calendar year.",
        )
        dna_durability = st.number_input(
            "DNA durability (years)", min_value=1, max_value=10_000,
            key="dna_durability",
            help="Years before the archive must be synthesized again. No replacement occurs at the exact end of the horizon.",
        )

    with st.container(key="cost_model_amazon"):
        _render_model_reset("amazon")
        st.caption("Price reference")
        amazon_base_year = st.number_input(
            "Amazon price base year", min_value=2000, max_value=2500,
            key="amazon_base_year",
            help="Calendar year to which all Amazon prices below apply.",
        )
        amazon_decline = st.number_input(
            "Amazon annual price decline (%)", min_value=0.0, max_value=99.99,
            key="amazon_decline",
            help="Annual reduction applied to Amazon request, retrieval, and storage prices.",
        )
        st.caption("Base-year prices")
        amazon_put_per_1000 = st.number_input(
            "Write requests (USD/1,000)", min_value=0.0,
            format="%.6f",
            key="amazon_put_per_1000",
            help="Charge for 1,000 requests when the archive is initially written.",
        )
        amazon_restore_per_1000 = st.number_input(
            "Bulk restore requests (USD/1,000)", min_value=0.0,
            format="%.6f",
            key="amazon_restore_per_1000",
            help="Charge for 1,000 bulk restore-job requests. Asset size determines the request count.",
        )
        amazon_retrieval_per_tb = st.number_input(
            "Bulk data retrieval (USD/TB)", min_value=0.0,
            format="%.6f",
            key="amazon_retrieval_per_tb",
            help="Capacity charge for retrieving one TB of archived data.",
        )
        amazon_storage_per_tb_month = st.number_input(
            "Storage (USD/TB/month)", min_value=0.0,
            format="%.6f",
            key="amazon_storage_per_tb_month",
            help="Recurring monthly charge to retain one TB in Deep Archive.",
        )

    with st.container(key="cost_model_azure"):
        _render_model_reset("azure")
        st.caption("Price reference")
        azure_base_year = st.number_input(
            "Azure price base year", min_value=2000, max_value=2500,
            key="azure_base_year",
            help="Calendar year to which all Azure prices below apply.",
        )
        azure_decline = st.number_input(
            "Azure annual price decline (%)", min_value=0.0, max_value=99.99,
            key="azure_decline",
            help="Annual reduction applied to Azure request, retrieval, and storage prices.",
        )
        st.caption("Base-year prices")
        azure_write_per_1000 = st.number_input(
            "Write requests (USD/1,000)", min_value=0.0,
            format="%.6f",
            key="azure_write_per_1000",
            help="Charge for 1,000 requests when the archive is initially written.",
        )
        azure_read_per_1000 = st.number_input(
            "Read requests (USD/1,000)", min_value=0.0,
            format="%.6f",
            key="azure_read_per_1000",
            help="Charge for 1,000 retrieval requests. Asset size determines the request count.",
        )
        azure_retrieval_per_tb = st.number_input(
            "Data retrieval (USD/TB)", min_value=0.0,
            format="%.6f",
            key="azure_retrieval_per_tb",
            help="Capacity charge for retrieving one TB from the Archive tier.",
        )
        azure_storage_per_tb_month = st.number_input(
            "Storage (USD/TB/month)", min_value=0.0,
            format="%.6f",
            key="azure_storage_per_tb_month",
            help="Recurring monthly charge to retain one TB in the Archive tier.",
        )

    with st.container(key="cost_model_tape"):
        _render_model_reset("tape")
        st.caption("Price reference")
        tape_base_year = st.number_input(
            "Tape price base year", min_value=2000, max_value=2500,
            key="tape_base_year",
            help="Calendar year to which the tape media, hardware, and energy prices apply.",
        )
        tape_durability = st.number_input(
            "Tape durability (years)", min_value=1, max_value=1_000,
            key="tape_durability",
            help="Years between complete tape media replacement writes.",
        )
        st.caption("Base-year prices")
        st.caption(
            "Tape cartridge and hardware values are added together. If your hardware estimate already "
            "includes cartridges/media, set tape cartridges to 0 to avoid double counting."
        )
        tape_media_per_tb = st.number_input(
            "Tape cartridges (USD/TB per write)", min_value=0.0,
            format="%.6f",
            key="tape_media_per_tb",
            help=(
                "Physical tape cartridge/media cost for one TB. The model charges it on the initial "
                "archive write and on each replacement write. Set this to 0 if the hardware value "
                "already includes cartridges."
            ),
        )
        tape_hardware_per_tb = st.number_input(
            "Tape library/drives (USD/TB amortized)", min_value=0.0,
            format="%.6f",
            key="tape_hardware_per_tb",
            help=(
                "Tape library, drive, robotics, and supporting hardware cost allocated per TB. The "
                "model spreads this over the selected tape durability period as annual maintenance."
            ),
        )
        tape_energy_per_tb_year = st.number_input(
            "Energy (USD/TB/year)", min_value=0.0,
            format="%.6f",
            key="tape_energy_per_tb_year",
            help="Annual energy cost to retain one TB in the tape system.",
        )
        st.caption("Annual price declines")
        tape_media_decline = st.number_input(
            "Tape cartridge decline (%)", min_value=0.0, max_value=99.99,
            key="tape_media_decline",
            help="Annual reduction applied to tape cartridge/media purchase prices.",
        )
        tape_hardware_decline = st.number_input(
            "Tape library/drives decline (%)", min_value=0.0, max_value=99.99,
            key="tape_hardware_decline",
            help="Annual reduction applied to amortized tape library, drive, and robotics costs.",
        )
        tape_energy_decline = st.number_input(
            "Tape energy decline (%)", min_value=0.0, max_value=99.99,
            key="tape_energy_decline",
            help="Annual reduction applied to tape energy costs.",
        )

    with st.container(key="cost_model_custom"):
        _render_model_reset("custom")
        custom_name = st.text_input(
            "Display name", key="custom_name",
            help="Name used for the custom technology in charts, tables, and downloads.",
        )
        custom_base_year = st.number_input(
            "Price base year", min_value=2000, max_value=2500,
            key="custom_base_year", help="Year to which all custom prices apply.",
        )
        custom_write_tb = st.number_input(
            "Initial write cost (USD/TB)", min_value=0.0,
            key="custom_write_tb", help="Capacity-based charge to write or replace one TB.",
        )
        custom_write_asset = st.number_input(
            "Write request cost (USD/asset)", min_value=0.0,
            key="custom_write_asset",
            help="Per-file or per-object charge applied when the archive is written or replaced.",
        )
        custom_storage_tb_year = st.number_input(
            "Annual storage cost (USD/TB)", min_value=0.0,
            key="custom_storage_tb_year",
            help="Recurring cost to retain one TB for one year.",
        )
        custom_retrieval_tb = st.number_input(
            "Retrieval cost (USD/TB)", min_value=0.0,
            key="custom_retrieval_tb",
            help="Capacity-based charge for each TB retrieved.",
        )
        custom_retrieval_asset = st.number_input(
            "Retrieval request cost (USD/asset)", min_value=0.0,
            key="custom_retrieval_asset",
            help="Per-file or per-object charge for the expected assets retrieved each year.",
        )
        custom_decline = st.number_input(
            "Annual price decline (%)", min_value=0.0, max_value=99.99,
            key="custom_decline",
            help="Annual percentage reduction applied to every custom price.",
        )
        custom_replacement = st.number_input(
            "Replacement interval (years)", min_value=0, max_value=10_000,
            key="custom_replacement",
            help="Years between complete rewrites. Use 0 for a service with no replacement writes.",
        )

    # Phones: the panel carries its own Calculate bar (sticky at the sheet's
    # bottom) so a cost change commits without closing the panel. Hidden on
    # desktop, where the header button covers it.
    calculate_panel = st.form_submit_button(
        "Calculate",
        key="calculate_panel",
        type="primary",
        width="stretch",
    )

_main_dg._form_data = None

submitted = calculate_header or calculate_scenario or calculate_panel

if submitted:
    try:
        submitted_widgets = _snapshot_widgets()
        candidate = _scenario_from_widgets(submitted_widgets)
    except (KeyError, ValueError) as error:
        st.error(str(error))
        st.stop()
    st.session_state["committed_widgets"] = submitted_widgets
    params = candidate.to_query_params()
    params.update(
        {
            "projection_end": str(int(submitted_widgets["projection_end"])),
            "log_scale": str(bool(submitted_widgets["log_scale"])),
            "theme": theme,
        }
    )
    st.query_params.from_dict(params)
    # Phones: the next rerun pops this flag and renders the iframe that
    # closes the sidebar sheet so the committed results are visible.
    st.session_state["close_sidebar_on_mobile"] = True
    st.rerun()

committed_widgets = st.session_state["committed_widgets"]

# Rendered only on the rerun after Calculate (see _SIDEBAR_CLOSE_JS).
if st.session_state.pop("close_sidebar_on_mobile", False):
    st.iframe(_SIDEBAR_CLOSE_JS, width=1, height=1)

# All graphs render from the last calculated inputs; widget edits live inside
# the form and do not rerun the script, so the graphs only change when a
# Calculate submit commits a new snapshot.
try:
    scenario = _scenario_from_widgets(committed_widgets)
except ValueError as error:
    st.error(str(error))
    st.stop()
if not scenario.technologies:
    st.error("Select at least one storage technology.")
    st.stop()
projection_end = int(committed_widgets["projection_end"])
log_scale = bool(committed_widgets["log_scale"])

result = _cached_simulation(scenario)
projection = _cached_projection(scenario, int(projection_end))
dna_curve_end = max(
    scenario.dna_cost_base_year,
    min(2500, scenario.start_year + scenario.horizon_years - 1),
)
dna_costs = _cached_dna_costs(scenario, dna_curve_end)
use_present_value = scenario.discount_rate_percent > 0
value_column = "present_value_usd" if use_present_value else "total_cost_usd"
chart_palette = palette_for(theme)

st.markdown('<div class="page-kicker">Archival storage economics</div>', unsafe_allow_html=True)
st.title("DNA Storage Cost Explorer")
st.markdown(
    '<p class="page-deck">Compare the long-run cost of DNA, cloud archive, tape, and custom storage '
    'under one consistent workload.</p>',
    unsafe_allow_html=True,
)
st.markdown(
    f"""
    <div class="model-strip">
        <span class="model-item"><strong>Model v{result.metadata['model_version']}</strong></span>
        <span class="model-item">{result.metadata['currency']}</span>
        <span class="model-item">Reviewed {result.metadata['last_reviewed']}</span>
        <span class="model-item">{result.metadata['disclaimer']}</span>
        <span class="model-item"><a href="https://doi.org/10.48550/arXiv.2608.26342" target="_blank" rel="noopener">Paper / About</a></span>
    </div>
    """,
    unsafe_allow_html=True,
)

totals = result.totals.sort_values(value_column)
cheapest = totals.iloc[0]
cheapest_name = str(cheapest["technology"])
cheapest_label = cheapest_name if len(cheapest_name) <= 26 else f"{cheapest_name[:25]}..."
dna_rows = totals[totals["technology"] == "DNA"]
dna_total = float(dna_rows.iloc[0][value_column]) if not dna_rows.empty else None
period_end = scenario.start_year + scenario.horizon_years - 1
cost_basis = "Present value" if use_present_value else "Undiscounted"
scale_label = "Log scale" if log_scale else "Linear scale"
st.markdown(
    f"""
    <div class="scenario-bar">
        <span class="scenario-label">Active scenario</span>
        <span><strong>{scenario.start_year}-{period_end}</strong> archive period</span>
        <span><strong>{len(scenario.technologies)}</strong> technologies</span>
        <span>{cost_basis}</span>
        <span>{scale_label}</span>
    </div>
    """,
    unsafe_allow_html=True,
)

metric_columns = st.columns(4)
metric_columns[0].metric("Archive", f"{scenario.archive_size_tb:,.3g} TB")
metric_columns[1].metric(
    "Data objects",
    _quantity(scenario.number_of_assets),
    help=f"{scenario.number_of_assets:,.0f} total data objects",
)
metric_columns[2].metric(
    f"Lowest: {cheapest_label}",
    _money(float(cheapest[value_column])),
    help=f"Lowest lifecycle cost: {cheapest_name}",
)
if dna_total is None:
    metric_columns[3].metric("Annual retrieval", f"{scenario.annual_retrieval_percent:,.2f}%")
else:
    lowest_cost = float(cheapest[value_column])
    if lowest_cost > 0:
        ratio = dna_total / lowest_cost
        metric_columns[3].metric(
            "DNA lifecycle cost", _money(dna_total), f"{ratio:,.2g}x lowest", delta_color="off"
        )
    else:
        metric_columns[3].metric("DNA lifecycle cost", _money(dna_total))

st.markdown('<div class="workspace-kicker">Analysis workspace</div>', unsafe_allow_html=True)
overview_tab, outlook_tab, dna_cost_tab, assumptions_tab, about_tab = st.tabs(
    ["Lifecycle", "Start-year outlook", "DNA unit costs", "Assumptions", "About"]
)

with overview_tab:
    _section_intro(
        "Lifecycle comparison",
        "Cumulative lifecycle cost for one archive opened in the selected start year. "
        "Each line includes initial and replacement writes, storage or operation, and expected retrieval.",
    )
    lifecycle_figure = lifecycle_chart(result, use_present_value, log_scale, theme=theme)
    _plotly_chart_with_placeholder(
        lifecycle_figure,
        key="chart_lifecycle",
        filename="dna-storage-lifecycle",
    )
    _chart_downloads(
        "lifecycle",
        result.yearly.to_csv(index=False).encode("utf-8"),
        "dna-storage-lifecycle",
        "chart_lifecycle",
    )
    st.markdown(
        """
        <div class="chart-divider">
            <h3>Cost composition</h3>
            <p>Undiscounted write and replacement, retrieval, and recurring storage or operating costs.</p>
        </div>
        """,
        unsafe_allow_html=True,
    )
    breakdown_figure = breakdown_chart(result, log_scale, theme=theme)
    _plotly_chart_with_placeholder(
        breakdown_figure,
        key="chart_breakdown",
        filename="dna-storage-cost-components",
    )
    breakdown_columns = ["technology", *[component for component in (
        "write_cost_usd", "read_cost_usd", "maintenance_cost_usd", "total_cost_usd"
    )]]
    _chart_downloads(
        "breakdown",
        result.totals[breakdown_columns].to_csv(index=False).encode("utf-8"),
        "dna-storage-cost-components",
        "chart_breakdown",
    )

with outlook_tab:
    _section_intro(
        "Start-year sensitivity",
        "The same archive workload and retention horizon are recalculated for every possible storage "
        "start year. A crossover is the first start year for which DNA's lifecycle cost is no greater "
        "than the comparison technology.",
    )
    projection_figure = projection_chart(projection, use_present_value, log_scale, theme=theme)
    _plotly_chart_with_placeholder(
        projection_figure,
        key="chart_projection",
        filename="dna-storage-start-year-outlook",
    )
    _chart_downloads(
        "projection",
        projection.to_csv(index=False).encode("utf-8"),
        "dna-storage-start-year-outlook",
        "chart_projection",
    )
    st.markdown(
        """
        <div class="chart-divider">
            <h3>Crossover years</h3>
            <p>First modeled start year in which DNA reaches or undercuts each comparison technology.</p>
        </div>
        """,
        unsafe_allow_html=True,
    )
    crossovers = find_crossover_years(projection)
    if crossovers:
        crossover_rows = [
            {
                "Comparison": f"DNA <= {technology}",
                "First start year": str(year) if year is not None else f"Not by {projection_end}",
            }
            for technology, year in crossovers.items()
        ]
        _themed_table(
            ["Comparison", "First start year"],
            [[row["Comparison"], row["First start year"]] for row in crossover_rows],
        )
    else:
        st.info("Include DNA and at least one comparison technology to calculate crossover years.")

with dna_cost_tab:
    _section_intro(
        "DNA unit economics",
        f"Both curves begin with the editable {scenario.dna_cost_base_year} unit costs and apply the "
        f"selected annual decline rates through {dna_curve_end}. They are unit-cost assumptions, "
        "not lifecycle totals.",
    )
    chart_columns = st.columns(2)
    with chart_columns[0]:
        synthesis_title = "DNA synthesis cost trajectory"
        synthesis_figure = dna_unit_cost_chart(
            dna_costs,
            "synthesis_cost_usd_per_mb",
            synthesis_title,
            chart_palette["unit_cost_colors"]["synthesis"],
            log_scale,
            theme=theme,
        )
        _plotly_chart_with_placeholder(
            synthesis_figure,
            key="chart_dna_synthesis",
            filename="dna-synthesis-cost-trajectory",
        )
        st.caption(
            "Modeled cost to synthesize enough DNA bases for 1 MB of logical data. "
            "Redundancy, indexing, and coding overhead are not added here."
        )
        synthesis_data = dna_costs[["year", "synthesis_cost_usd_per_mb"]]
        _chart_downloads(
            "dna_synthesis",
            synthesis_data.to_csv(index=False).encode("utf-8"),
            "dna-synthesis-cost-trajectory",
            "chart_dna_synthesis",
        )
    with chart_columns[1]:
        sequencing_title = "DNA sequencing cost trajectory"
        sequencing_figure = dna_unit_cost_chart(
            dna_costs,
            "sequencing_cost_usd_per_mb",
            sequencing_title,
            chart_palette["unit_cost_colors"]["sequencing"],
            log_scale,
            theme=theme,
        )
        _plotly_chart_with_placeholder(
            sequencing_figure,
            key="chart_dna_sequencing",
            filename="dna-sequencing-cost-trajectory",
        )
        st.caption(
            "Modeled cost to sequence and retrieve 1 MB of logical data. "
            "It is applied to the share of the archive retrieved each year."
        )
        sequencing_data = dna_costs[["year", "sequencing_cost_usd_per_mb"]]
        _chart_downloads(
            "dna_sequencing",
            sequencing_data.to_csv(index=False).encode("utf-8"),
            "dna-sequencing-cost-trajectory",
            "chart_dna_sequencing",
        )

with assumptions_tab:
    assumptions = load_assumptions()
    _section_intro(
        "Assumptions and scope",
        "A concise record of the active scenario, included cost categories, and source references.",
    )
    contract = pd.DataFrame(
        [
            ("Archive period", f"{scenario.start_year}-{scenario.start_year + scenario.horizon_years - 1}"),
            ("Annual retrieval", f"{scenario.annual_assets_retrieved:,.2f} assets ({scenario.annual_retrieval_percent:,.2f}%)"),
            ("Discounting", f"{scenario.discount_rate_percent:,.2f}% real; present value at storage start"),
            ("Currency", result.metadata["currency"]),
            (
                "DNA synthesis baseline",
                f"${scenario.dna_synthesis_cost_per_mb:,.6g}/MB in {scenario.dna_cost_base_year}",
            ),
            (
                "DNA sequencing baseline",
                f"${scenario.dna_sequencing_cost_per_mb:,.6g}/MB in {scenario.dna_cost_base_year}",
            ),
            (
                "Amazon baseline",
                f"${scenario.amazon_storage_usd_per_mb_month * 1_000_000:,.6g}/TB/month storage; "
                f"${scenario.amazon_bulk_retrieval_usd_per_mb * 1_000_000:,.6g}/TB retrieval "
                f"in {scenario.amazon_price_base_year}",
            ),
            (
                "Azure baseline",
                f"${scenario.azure_storage_usd_per_mb_month * 1_000_000:,.6g}/TB/month storage; "
                f"${scenario.azure_retrieval_usd_per_mb * 1_000_000:,.6g}/TB retrieval "
                f"in {scenario.azure_price_base_year}",
            ),
            (
                "Tape baseline",
                f"${scenario.tape_media_usd_per_tb:,.6g}/TB media; "
                f"${scenario.tape_hardware_usd_per_tb:,.6g}/TB hardware; "
                f"${scenario.tape_energy_usd_per_tb_year:,.6g}/TB/year energy "
                f"in {scenario.tape_price_base_year}",
            ),
            ("Horizon convention", result.metadata["horizon_convention"]),
        ],
        columns=["Item", "Value"],
    )
    contract_column, context_column = st.columns([1.45, 1], gap="large")
    with contract_column:
        st.subheader("Scenario contract")
        _themed_table(
            list(contract.columns),
            [list(row) for row in contract.itertuples(index=False)],
        )
    with context_column:
        st.subheader("Cost scope")
        st.write(
            "Included: DNA synthesis and sequencing, cloud write/retrieval/storage charges, "
            "tape media/hardware/energy assumptions, and the selected custom capacity/request charges. "
            "Excluded: labor, cloud egress, taxes, retrieval latency, minimum-storage penalties, "
            "facilities, and unmodeled migrations."
        )
        st.subheader("Sources")
        for source in assumptions["sources"].values():
            st.markdown(f"- [{source['label']}]({source['url']})")

with about_tab:
    _section_intro(
        "About this explorer",
        "An interactive implementation of a DNA storage cost model for comparing long-run archival "
        "economics across DNA, cloud archive, tape, and user-defined storage systems.",
    )
    about_column, paper_column = st.columns([1.25, 1], gap="large")
    with about_column:
        st.subheader("Purpose")
        st.write(
            "The explorer turns the model assumptions into adjustable controls, then recalculates "
            "lifecycle cost, start-year sensitivity, and DNA unit-cost trajectories for the selected "
            "archive workload."
        )
        st.subheader("Model boundary")
        st.write(
            "The results are scenario estimates, not procurement quotes. They use the visible inputs "
            "and listed assumptions, and exclude operational details such as labor, taxes, cloud egress, "
            "retrieval latency, and migration execution risk."
        )
    with paper_column:
        st.subheader("Reference paper")
        st.markdown(
            "[DNA Storage Cost Model](https://doi.org/10.48550/arXiv.2608.26342)"
        )
        st.write(
            "Use the linked paper as the source reference for the model framing and baseline assumptions."
        )
        st.subheader("Contact")
        st.markdown("[alex@el-shaikh.com](mailto:alex@el-shaikh.com)")
