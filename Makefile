# ============================================================
# Makefile for generate_polychora
# Requires: GNU Make, .NET 8 SDK
# Usage:
#   make                  # build + generate topology for missing files
#   make vertices         # regenerate all vertex JSONs
#   make topology         # compute topology for all missing files
#   make topology/pen     # (re)compute topology for one polytope
#   make rebuild-topology # delete all topology files and regenerate
#   make check            # Euler-characteristic check on all topology files
#   make clean            # remove topology_output/ and build artefacts
# ============================================================

DOTNET       := dotnet
PROJ         := generate_polychora.csproj
VDIR         := vertex_output
TDIR         := topology_output

# ── auto-discover files ──────────────────────────────────────
VERTEX_JSONS := $(wildcard $(VDIR)/*.json)
TOPO_JSONS   := $(patsubst $(VDIR)/%.json, $(TDIR)/%.json, $(VERTEX_JSONS))
NAMES        := $(patsubst $(VDIR)/%.json, %, $(VERTEX_JSONS))

# ── phony targets ────────────────────────────────────────────
.PHONY: all vertices topology rebuild-topology check clean \
        $(addprefix topology/, $(NAMES))

# ── default ──────────────────────────────────────────────────
all: build topology

# ── build .NET project ───────────────────────────────────────
build: $(PROJ)
	$(DOTNET) build -q

# ── vertex generation ────────────────────────────────────────
vertices: build
	$(DOTNET) run -- vertices

# ── topology: generate all missing files ─────────────────────
topology: build $(TDIR) $(TOPO_JSONS)

# Pattern rule: each topology file depends on its vertex file.
# If the topology file is missing (or older), invoke topology-one.
$(TDIR)/%.json: $(VDIR)/%.json | $(TDIR)
	@echo "  → $*"
	$(DOTNET) run -- topology-one $*

$(TDIR):
	mkdir -p $(TDIR)

# ── per-name convenience: make topology/pen ──────────────────
$(addprefix topology/, $(NAMES)): topology/%:
	rm -f $(TDIR)/$*.json
	$(DOTNET) run -- topology-one $*

# ── rebuild everything ────────────────────────────────────────
rebuild-topology: | $(TDIR)
	rm -f $(TDIR)/*.json
	$(MAKE) topology

# ── integrity / Euler check ──────────────────────────────────
check:
	@python check_topology.py $(TDIR)

# ── clean ─────────────────────────────────────────────────────
clean:
	rm -rf $(TDIR) bin/ obj/
	$(DOTNET) clean -q

clean-topology:
	rm -f $(TDIR)/*.json

# ── status overview ───────────────────────────────────────────
status:
	@echo "vertex_output:   $$(ls $(VDIR)/*.json 2>/dev/null | wc -l) files"
	@echo "topology_output: $$(ls $(TDIR)/*.json 2>/dev/null | wc -l) files"
	@python -c "\
import os; \
v=set(os.listdir('$(VDIR)')); \
t=set(os.listdir('$(TDIR)') if os.path.isdir('$(TDIR)') else []); \
m=sorted(v-t); \
print('Missing topology: '+', '.join(f[:-5] for f in m) if m else 'All topology files present.')"
