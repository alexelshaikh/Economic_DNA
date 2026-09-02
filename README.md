# An Economic Analysis of DNA-based Data Storage Systems  
## Alex El-Shaikh, Bernhard Seeger, and Thomas Heinis

## Interactive cost explorer

The repository includes a journalist-facing web calculator in `streamlit_app.py`. It uses a
validated, numeric implementation of the paper model and exposes archive size, average asset
size, annual retrieval, start year, retention horizon, discounting, decline assumptions, and
technology selection. DNA synthesis and sequencing base-year costs are editable, and a custom
storage service can combine per-TB, per-asset, retrieval, recurring, decline, and replacement
assumptions. Results include lifecycle, start-year, synthesis, and sequencing charts, crossover
years, source metadata, and CSV/PNG/SVG downloads.

Run it locally:

```bash
pip install -r requirements.txt
streamlit run streamlit_app.py
```

The paper baseline is 1 TB stored as 1,000 assets of 1 GB each, with 1% of the assets retrieved
per year for 100 years. A 100-year horizon includes the start year through start year + 99.

### Docker Compose

Start the app locally from the repository directory:

```bash
docker compose up -d --build
```

Open `http://localhost:8501`. Later starts can omit `--build` unless the source or dependencies
changed. View status and logs, restart the app, or stop the stack with:

```bash
docker compose ps
docker compose logs -f app
docker compose restart app
docker compose down
```

Set `APP_PORT` in `.env` to use another local port.

### VPS deployment

Copy `.env.example` to `.env`, replace `dna-cost.example.org` with a domain whose DNS points to
the VPS, and run:

```bash
docker compose --profile production up -d --build
```

Caddy proxies the app and manages HTTPS automatically. The VPS must allow inbound TCP 80/443
and UDP 443. Keep `.env` out of version control. For embargoed access, configure Streamlit OIDC
or an authentication gateway before making the DNS record public.

Run all automated checks with:

```bash
python -m unittest discover -v
```

This repository contains the code and data accompanying the paper **“An Economic Analysis of DNA-based Data Storage Systems.”** It provides fully reproducible notebooks for all main-text and supplementary figures, along with scripts and utilities to fetch and cache datasets.

> **TL;DR:**  
> 1️⃣ Ensure you have an **internet connection**.  
> 2️⃣ `pip install -r requirements.txt`  
> 3️⃣ Open `Figures.ipynb` (or `Figures_Supplementary.ipynb`) and **run the very first cell once**.  
> 4️⃣ After that, any figure-generating cell can be run **independently**.

---

## 📂 Repository Structure

- `Figures.ipynb` — Reproduces all figures in the **main manuscript**. Markdown cells label each figure for easy navigation.  
- `Figures_Supplementary.ipynb` — Reproduces all **supplementary figures** referenced by the main text.  
- `requirements.txt` — Python dependencies for a clean, reproducible environment.  
- `data/` — Contains locally cached datasets. Some data are fetched on first run and stored here.  
- `models/` — Contains the storage cost models introduced in the paper.  
- `figs/` — Output directory for generated plots when using the helper function `save(plot_name)`.  
- `preamble.py` — Loads all required modules and packages. This is automatically imported by running the first cell in `Figures.ipynb` or `Figures_Supplementary.ipynb`.  
- `storage_service.py` — Provides methods to load specific storage models (e.g., tape, DNA).

---

## 🧰 Requirements

You will need **Python 3.10+** and the packages listed in `requirements.txt`. Install them with:

```bash
pip install -r requirements.txt
```

**Included dependencies:**

```
requests
sympy
pandas
seaborn
numpy
matplotlib
scipy
xlrd>=2.0.1
tqdm
jupyter-server
```

> ⚠️ **Internet required:** Some datasets and storage models are downloaded on first use. If you’re offline, these steps will fail.

---

## 🚀 Quick Start

1. **Clone** the repository:
   ```bash
   git clone https://github.com/alexelshaikh/Economic_DNA
   cd Economic_DNA
   ```

2. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

3. **Launch Jupyter**:
   ```bash
   jupyter lab
   # or
   jupyter notebook
   ```

4. **Open a notebook**:
   - Main figures: `Figures.ipynb`  
   - Supplementary figures: `Figures_Supplementary.ipynb`

5. **Run the first cell once** to set up paths, configuration, and download any required data.

6. **Generate figures**:  
   After initialization, any figure cell can be executed independently and in any order.

---

## 📊 Reproducing the Figures

- **Main manuscript**:  
  Open `Figures.ipynb`. Each figure is labeled via Markdown headings (e.g., “Figure 1”, “Figure 2”, …). Run the first cell, then execute the desired figure cells.

- **Supplementary figures**:  
  Open `Figures_Supplementary.ipynb`. Run the first cell, then jump to any supplementary figure section and execute.

> 💡 If you change paths or environment variables, re-run the first cell to refresh the session state.

---

## 💾 Data Access & Caching

- Datasets are **downloaded on demand** (via `requests`) the first time a relevant cell runs.  
- Downloaded files are **cached locally** (typically under `data/`) for faster, offline-friendly subsequent runs.

---

## 🧭 Troubleshooting

- **No internet / download errors:** Check your connection, then re-run the first cell.  
- **Excel reading errors:** Ensure `xlrd>=2.0.1` is installed.  
- **Jupyter not found:** Confirm `jupyter-server` is installed and your environment is activated.  
- **Permission issues in `data/` or `figs/`:** Make sure you have write permissions to the repo directory.

---

## 📝 Citation

If you use this code or data, please cite our paper.  
A BibTeX entry will be provided here once the paper is published.

---

## 📄 License

This project is licensed under the **MIT License** — see the [LICENSE](./LICENSE) file for details.

---

## 📬 Contact

For questions, issues, or requests, please open a GitHub issue or contact **Alex El-Shaikh** at:  
📧 `a.elshaikh@imperial.ac.uk`
