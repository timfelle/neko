ACLOCAL_AMFLAGS = -I m4 --install
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = neko.pc

SUBDIRS = src

if ENABLE_CONTRIB
SUBDIRS += contrib
endif

if ENABLE_PFUNIT
SUBDIRS += tests/unit
endif

if ENABLE_DOXYGEN
SUBDIRS += doc
endif

CLEANFILES = *.mod\
	     *.smod


bin_SCRIPTS = makeneko

man_MANS = \
	doc/man/neko.1\
	doc/man/makeneko.1

EXTRA_DIST = \
	makeneko.in\
	regen.sh\
	flinter_rc.yml\
	tests/bench/ax/setup.f90\
	tests/bench/ax/Makefile\
	tests/bench/ax/axkernel.f90\
	tests/bench/ax/driver.f90\
	tests/bench/gs/Makefile\
	tests/bench/gs/driver.f90\
	tests/bench/tgv32/tgv.case\
	tests/bench/tgv32/README.md\
	tests/bench/tgv32/extract.sh\
	tests/bench/tgv32/tgv.f90\
	examples/turb_channel/turb_channel.case\
	examples/turb_channel/README.md\
	examples/turb_channel/turb_channel.f90\
	examples/turb_channel/box.nmsh\
	examples/turb_pipe/turb_pipe.case\
	examples/turb_pipe/turb_pipe_source.case\
	examples/turb_pipe/turb_pipe.f90\
	examples/turb_pipe/turb_pipe.nmsh\
	examples/turb_pipe/turb_pipe.geo\
	examples/turb_pipe/README.md\
	examples/hemi/hemi.nmsh\
	examples/hemi/README.md\
	examples/hemi/hemi.case\
	examples/lid/lid.f90\
	examples/lid/lid2d.nmsh\
	examples/lid/lid.nmsh\
	examples/lid/README.md\
	examples/lid/lid2d.case\
	examples/lid/lid.case\
	examples/lid/elem_dist.csv\
	examples/lid/plot_stat.py\
	examples/poisson/ax_poisson.f90\
	examples/poisson/setup.f90\
	examples/poisson/Makefile\
	examples/poisson/README.md\
	examples/poisson/data/4096.nmsh\
	examples/poisson/data/16384.nmsh\
	examples/poisson/data/256.nmsh\
	examples/poisson/data/2048.nmsh\
	examples/poisson/data/128.nmsh\
	examples/poisson/data/32768.nmsh\
	examples/poisson/data/8192.nmsh\
	examples/poisson/data/512.nmsh\
	examples/poisson/data/1024.nmsh\
	examples/poisson/driver.f90\
	examples/2d_cylinder/2d_cylinder.nmsh\
	examples/2d_cylinder/2d_cylinder.case\
	examples/2d_cylinder/2d_cylinder.f90\
	examples/2d_cylinder/2d_cylinder.geo\
	examples/2d_cylinder/README.md\
	examples/cyl_boundary_layer/cyl_bl_rot.case\
	examples/cyl_boundary_layer/cyl_bl.f90\
	examples/cyl_boundary_layer/README.md\
	examples/cyl_boundary_layer/cyl_bl_user.case\
	examples/cyl_boundary_layer/cyl_bl_basic.case\
	examples/cyl_boundary_layer/cyl.nmsh\
	examples/rayleigh_benard/rayleigh.f90\
	examples/rayleigh_benard/README.md\
	examples/rayleigh_benard/box.nmsh\
	examples/rayleigh_benard/rayleigh.case\
	examples/tgv/262144.nmsh\
	examples/tgv/tgv.case\
	examples/tgv/tgv_expression.case\
	examples/tgv/README.md\
	examples/tgv/32768.nmsh\
	examples/tgv/tgv.f90\
	examples/tgv/512.nmsh\
	examples/rayleigh_benard_cylinder/probes.csv\
	examples/rayleigh_benard_cylinder/rayleigh.f90\
	examples/rayleigh_benard_cylinder/README.md\
	examples/rayleigh_benard_cylinder/rayleigh.nmsh\
	examples/rayleigh_benard_cylinder/rayleigh.geo\
	examples/rayleigh_benard_cylinder/rayleigh.case\
	examples/rayleigh_benard_cylinder/postprocess_probes.py\
	examples/rebounding_particles/README.md\
	examples/rebounding_particles/particles_init.py\
	examples/rebounding_particles/rebounding_particles_csv.case\
	examples/rebounding_particles/rebounding_particles_h5.case\
	examples/rebounding_particles/run.sh\
	examples/rebounding_particles/trajectory_plot_csv.py\
	examples/rebounding_particles/trajectory_plot_h5.py\
	examples/api/pyneko/README.md\
	examples/api/pyneko/cylinder.case\
	examples/api/pyneko/cylinder.py\
	examples/api/pyneko/cylinder.nmsh\
	examples/api/c/README.md\
	examples/api/c/cylinder.case\
	examples/api/c/cylinder.c\
	examples/api/c/cylinder.nmsh\
	tests/reframe/checks.py\
	tests/reframe/settings.py\
	tests/reframe/README.md\
	tests/reframe/src/get-tgv-dns.sh\
	tests/reframe/src/minihemi.case.template\
	tests/reframe/src/rayleigh.case.template\
	tests/reframe/src/rayleigh.f90\
	tests/reframe/src/rocm_select_gpu_device.sh\
	tests/reframe/src/tgv.case.template\
	tests/reframe/src/tgv.f90\
	tests/regression/valgrind/README.md\
	tests/regression/valgrind/run_valgrind_check.sh\
	tests/regression/valgrind/check_valgrind.py\
	tests/regression/valgrind/valgrind_runtime.supp\
	tests/integration/default_cpu_launcher.sh\
	tests/integration/default_cuda_launcher.sh\
	tests/integration/conftest.py\
	tests/integration/testlib.py\
	tests/integration/pytest.ini\
	tests/integration/meshes/small_test_cyl.nmsh\
	tests/integration/tests/test_demo/test_suite.py\
	tests/integration/tests/test_examples/test_examples.py\
	tests/integration/tests/test_cylinder/test_cylinder.py\
	tests/integration/tests/test_cylinder/cylinder_part1.case\
	tests/integration/tests/test_cylinder/reflog1_cpu_dp.csv\
	tests/integration/tests/test_cylinder/reflog1_cpu_sp.csv\
	tests/integration/tests/test_cylinder/reflog2_cpu_dp.csv\
	tests/integration/tests/test_cylinder/reflog2_cpu_sp.csv\
	tests/integration/tests/test_cylinder/reflog1_cuda_dp.csv\
	tests/integration/tests/test_cylinder/reflog1_cuda_sp.csv\
	tests/integration/tests/test_cylinder/reflog2_cuda_dp.csv\
	tests/integration/tests/test_cylinder/reflog2_cuda_sp.csv\
	tests/integration/tests/test_cylinder/reflog1_hip_dp.csv\
	tests/integration/tests/test_cylinder/reflog1_hip_sp.csv\
	tests/integration/tests/test_cylinder/reflog2_hip_dp.csv\
	tests/integration/tests/test_cylinder/reflog2_hip_sp.csv\
	tests/integration/tests/test_cylinder/reflog1_opencl_dp.csv\
	tests/integration/tests/test_cylinder/reflog1_opencl_sp.csv\
	tests/integration/tests/test_cylinder/reflog2_opencl_dp.csv\
	tests/integration/tests/test_cylinder/reflog2_opencl_sp.csv\
	patches/nvhpc_bge.patch\
	patches/json-fortran71.patch\
	patches/cce_stack.patch\
	patches/cce_time_state.patch\
	patches/fujitsu_memleak.patch

if ENABLE_FLINT
# Score every Fortran source individually instead of handing the whole src/
# tree to a single flint invocation. flint's parser allocates without bound on
# some sources, so the whole-tree call is reliably killed by the OOM killer,
# which in turn kills `make check`.
#
# This is only loosely comparable to .github/workflows/check_lint.yml: that
# workflow's per-file loop runs over the *changed* files of a pull request,
# not the whole tree, and hardcodes a skip list rather than a memory limit.
# The cost of covering every file here is real -- roughly two minutes and one
# line of output per source -- which is why nothing is gated on the result.
#
# Reporting only: like the plain `flint score` it replaces, this target never
# fails the build on code quality, whatever the scores are. 174 of the 482
# sources under src/ are below 10 today, and gating on that is a separate
# decision from fixing the crash.
#
# Sources generated at configure time are skipped, by ignoring any file that
# has a matching *.in template, so the set of files scored does not depend on
# whether this is an in-tree or a VPATH build.
#
# Both knobs below are deliberately plain make variables rather than configure
# substitutions: they are developer conveniences, overridable per invocation
# with e.g. `make lint FLINT_SCORE_GOAL=9`, and nothing in the build depends
# on their values.
#
# FLINT_SCORE_GOAL  Per-file score at or above which a file is considered
#                   fine; anything below is listed as improvable at the end of
#                   the run. Matches the bar check_lint.yml reports against.
# FLINT_MEM_LIMIT   Per-file address-space limit, in KiB. flint cannot score a
#                   handful of sources at all (it grows past 20 GiB and still
#                   does not finish); the limit turns that from an OOM kill of
#                   the whole build into a reported, non-fatal skip. Platforms
#                   that do not enforce RLIMIT_AS warn once and carry on.
FLINT_SCORE_GOAL = 10
FLINT_MEM_LIMIT = 2097152

# Sources flint cannot score at any memory limit -- it grows past 20 GiB on
# these and still does not terminate. Skipped outright so the run does not
# spend minutes reaching the limit; check_lint.yml skips the same two.
FLINT_EXCLUDE = src/adt/stack.f90 src/adt/htable.f90

.PHONY: lint
lint:
	@tmpdir=`mktemp -d 2> /dev/null || echo ./lint-tmp.$$$$`; \
	mkdir -p "$$tmpdir" || exit 1; \
	trap 'rm -rf "$$tmpdir"' EXIT HUP INT TERM; \
	lim='$(FLINT_MEM_LIMIT)'; \
	if ( ulimit -v $$lim ) 2> /dev/null; then :; else \
	  printf '%s\n' "Warning: this platform will not honour a $$lim KiB"; \
	  printf '%s\n' "address-space limit, so a pathological source can" \
	    "still exhaust memory and abort the build."; \
	  lim=unlimited; \
	fi; \
	find "$(srcdir)/src" \( -name '*.f90' -o -name '*.F90' \) \
	  | sort > "$$tmpdir/files"; \
	: > "$$tmpdir/improvable"; \
	: > "$$tmpdir/skipped"; \
	while IFS= read -r f; do \
	  if test -f "$$f.in"; then continue; fi; \
	  excluded=no; \
	  for x in $(FLINT_EXCLUDE); do \
	    case "$$f" in */$$x) excluded=yes;; esac; \
	  done; \
	  if test "$$excluded" = yes; then \
	    printf '  - %s: skipped, known to exhaust memory\n' "$$f"; \
	    printf '  - %s: on FLINT_EXCLUDE\n' "$$f" >> "$$tmpdir/skipped"; \
	    continue; \
	  fi; \
	  out=`( ulimit -v $$lim 2> /dev/null; \
	         $(FLINT) score -d 0 -r "$(srcdir)/flinter_rc.yml" "$$f" 2>&1 \
	       )`; \
	  status=$$?; \
	  score=`printf '%s\n' "$$out" | sed -n 's/.*-->|\([^|]*\)|<--.*/\1/p'`; \
	  if test -z "$$score"; then \
	    reason=`printf '%s\n' "$$out" | tail -n 1`; \
	    printf '  - %s: not scored, flint exited %s\n' "$$f" "$$status"; \
	    printf '  - %s: exit %s, %s\n' "$$f" "$$status" "$$reason" \
	      >> "$$tmpdir/skipped"; \
	  else \
	    printf '  - %s: %s\n' "$$f" "$$score"; \
	    if awk "BEGIN { exit !($$score < $(FLINT_SCORE_GOAL)) }"; then \
	      printf '  - %s: %s\n' "$$f" "$$score" >> "$$tmpdir/improvable"; \
	    fi; \
	  fi; \
	done < "$$tmpdir/files"; \
	if test -s "$$tmpdir/skipped"; then \
	  printf '%s\n' \
	    "Files flint could not score (see FLINT_EXCLUDE, FLINT_MEM_LIMIT):"; \
	  cat "$$tmpdir/skipped"; \
	fi; \
	if test -s "$$tmpdir/improvable"; then \
	  printf 'Files that can be improved (score below %s):\n' \
	    '$(FLINT_SCORE_GOAL)'; \
	  cat "$$tmpdir/improvable"; \
	  printf 'Run "%s lint -r %s/flinter_rc.yml <file>" for the\n' \
	    '$(FLINT)' '$(srcdir)'; \
	  printf 'suggested improvements.\n'; \
	fi; \
	exit 0

check-local: lint

else
.PHONY: lint
lint:
	@echo "flint support is disabled: nothing to lint."
endif

if ENABLE_PFUNIT
.PHONY: refresh-unit-test-makefiles
refresh-unit-test-makefiles:
	$(MAKE) -C tests/unit refresh-unit-test-makefiles-local
else
.PHONY: refresh-unit-test-makefiles
refresh-unit-test-makefiles:
	@echo "pFUnit support is disabled: no unit-test Makefiles to refresh."
endif
