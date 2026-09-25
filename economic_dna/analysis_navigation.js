(() => {
  const win = window;
  const doc = document;
  win.__dnaAnalysisNavigation?.dispose();
  let pending = null;
  let frame = 0;
  let sequence = 0;

  const finish = () => {
    win.cancelAnimationFrame(frame);
    frame = 0;
    if (!pending) return;
    win.clearTimeout(pending.noticeTimer);
    pending.preview.remove();
    pending.notice.remove();
    pending.root.classList.remove("analysis-transition");
    pending.root.style.removeProperty("min-height");
    pending.root.removeAttribute("aria-busy");
    pending = null;
  };

  // One inert visual snapshot bridges the server/Plotly render. It contains
  // no live widgets, export handlers, iframes, or duplicate element IDs.
  const snapshot = panel => {
    const copy = panel.cloneNode(true);
    copy.querySelectorAll('[data-testid="stElementContainer"]:has(iframe), [data-testid="stElementContainer"]:has(.analysis-ready)')
      .forEach(node => node.remove());
    copy.querySelectorAll("iframe, script, style, .analysis-ready").forEach(node => node.remove());
    const columns = panel.querySelectorAll('[data-testid="stColumn"]');
    copy.querySelectorAll('[data-testid="stColumn"]').forEach((column, index) => {
      const width = `${columns[index].getBoundingClientRect().width}px`;
      Object.assign(column.style, {flex: `0 0 ${width}`, minWidth: width, maxWidth: width});
    });
    const buttons = panel.querySelectorAll("button");
    copy.querySelectorAll("button").forEach((button, index) => {
      const style = win.getComputedStyle(buttons[index]);
      for (const property of ["font-size", "color", "height", "min-height"]) {
        button.style.setProperty(property, style.getPropertyValue(property));
      }
    });
    const nodes = [copy, ...copy.querySelectorAll("*")];
    const ids = new Map();
    const prefix = `analysis-preview-${++sequence}-`;
    for (const node of nodes) {
      if (node.id) ids.set(node.id, prefix + node.id);
    }
    for (const node of nodes) {
      if (node.id) node.id = ids.get(node.id);
      for (const name of [...node.classList]) {
        if (name.startsWith("st-key-")) node.classList.remove(name);
      }
      for (const attribute of [...node.attributes]) {
        if (attribute.name.startsWith("on") || ["role", "tabindex", "name", "form"].includes(attribute.name)) {
          node.removeAttribute(attribute.name);
        } else {
          let value = attribute.value.replace(/url\(#([^)]+)\)/g, (match, id) => ids.has(id) ? `url(#${ids.get(id)})` : match);
          if (["href", "xlink:href"].includes(attribute.name) && value.startsWith("#") && ids.has(value.slice(1))) {
            value = "#" + ids.get(value.slice(1));
          }
          if (value !== attribute.value) node.setAttribute(attribute.name, value);
        }
      }
    }
    copy.removeAttribute("data-testid");
    copy.classList.add("analysis-preview");
    copy.setAttribute("aria-hidden", "true");
    copy.inert = true;
    return copy;
  };

  const ready = () => {
    frame = 0;
    if (!pending) return;
    const {root, target, generation, started} = pending;
    const selected = root.querySelector('[role="tab"][aria-selected="true"]');
    const panel = root.querySelector('[role="tabpanel"]');
    const marker = panel?.querySelector(".analysis-ready");
    const charts = panel ? [...panel.querySelectorAll('[class*="st-key-chart_"]')] : [];
    const complete = selected?.dataset.key === target
      && Number(marker?.dataset.renderId) > generation
      && !panel.querySelector('[data-stale="true"]')
      && charts.every(element => {
        const plot = element.querySelector(".js-plotly-plot");
        return plot?._fullLayout?.width > 0 && plot.querySelector(".main-svg");
      });
    pending.readyFrames = complete ? pending.readyFrames + 1 : 0;
    if (pending.readyFrames >= 2 || !root.isConnected || performance.now() - started > 30000) {
      finish();
    } else {
      frame = win.requestAnimationFrame(ready);
    }
  };

  const begin = tab => {
    const root = tab?.closest(".st-key-analysis_tabs");
    if (!root || tab.getAttribute("aria-selected") === "true") return;
    const target = tab.dataset.key;
    if (pending) {
      pending.target = target;
      pending.notice.textContent = `Loading ${tab.textContent.trim()}...`;
      pending.readyFrames = 0;
      return;
    }
    const panel = root.querySelector('[role="tabpanel"]');
    const marker = panel?.querySelector(".analysis-ready");
    if (!marker) return;
    const bounds = root.getBoundingClientRect();
    const top = panel.getBoundingClientRect().top - bounds.top;
    const preview = snapshot(panel);
    const notice = doc.createElement("div");
    notice.className = "analysis-loading-status";
    notice.setAttribute("role", "status");
    notice.textContent = `Loading ${tab.textContent.trim()}...`;
    notice.hidden = true;
    preview.style.top = `${top}px`;
    notice.style.top = `${top - 6}px`;
    root.style.minHeight = `${bounds.height}px`;
    root.append(preview, notice);
    root.classList.add("analysis-transition");
    root.setAttribute("aria-busy", "true");
    pending = {
      root, preview, notice, target, generation: Number(marker.dataset.renderId),
      started: performance.now(), readyFrames: 0,
      noticeTimer: win.setTimeout(() => { notice.hidden = false; }, 150),
    };
    frame = win.requestAnimationFrame(ready);
  };

  const click = event => begin(event.target.closest('[role="tab"]'));
  const pointerdown = event => {
    if (event.button === 0) begin(event.target.closest('[role="tab"]'));
  };
  const pointerup = event => {
    const selected = pending?.root.querySelector('[role="tab"][aria-selected="true"]');
    if (pending && selected?.dataset.key !== pending.target && !event.target.closest('.st-key-analysis_tabs [role="tab"]')) finish();
  };
  const keydown = event => {
    if (event.key === "Escape") finish();
    const tab = event.target.closest('.st-key-analysis_tabs [role="tab"]');
    if (!tab || !["ArrowLeft", "ArrowRight", "Home", "End"].includes(event.key)) return;
    const tabs = [...tab.closest('[role="tablist"]').querySelectorAll('[role="tab"]')];
    let index = tabs.indexOf(tab);
    if (event.key === "Home") index = 0;
    else if (event.key === "End") index = tabs.length - 1;
    else index = (index + (event.key === "ArrowRight" ? 1 : -1) + tabs.length) % tabs.length;
    begin(tabs[index]);
  };
  doc.addEventListener("click", click, true);
  doc.addEventListener("pointerdown", pointerdown, true);
  doc.addEventListener("pointerup", pointerup, true);
  doc.addEventListener("pointercancel", finish, true);
  doc.addEventListener("keydown", keydown, true);
  win.addEventListener("resize", finish);
  win.__dnaAnalysisNavigation = {
    dispose() {
      finish();
      doc.removeEventListener("click", click, true);
      doc.removeEventListener("pointerdown", pointerdown, true);
      doc.removeEventListener("pointerup", pointerup, true);
      doc.removeEventListener("pointercancel", finish, true);
      doc.removeEventListener("keydown", keydown, true);
      win.removeEventListener("resize", finish);
    },
  };
})();
