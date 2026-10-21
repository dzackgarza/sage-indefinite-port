# The Sage the acceptance workflow runs in: the pinned sagemath-dev build,
# rebuilt on Python 3.14 (dzack-research requires >=3.14) exactly as
# ~/ai-review-ci/scripts/sage-docker rebuilds the local sage-dev-develop
# container. Built and pushed by .github/workflows/sage-image.yml.
FROM sagemath/sagemath-dev@sha256:2928e82688bf4f182ae1f9a03187cfdb2b4ea4c4f0e902963c8130c1959666d9 AS build
COPY --from=ghcr.io/astral-sh/uv:0.11.17 /uv /usr/local/bin/uv
USER sage
WORKDIR /home/sage/sage
RUN set -eux; \
    uv python install 3.14.5; \
    python314="$(uv python find --managed-python --resolve-links 3.14.5)"; \
    uv pip install --python "$python314" --break-system-packages setuptools; \
    ./configure --with-python="$python314"; \
    make -j"$(nproc)" build; \
    test "$(./sage -python -c 'import sys; print(sys.version_info[:2])')" = "(3, 14)"; \
    rm -rf upstream; \
    uv cache clean; \
    MESONPY_EDITABLE_VERBOSE=1 ./sage -python -c 'import sage.all'

# One layer, so the source tarballs removed above do not stay in the image
# underneath. The Python 3.12 venv stays: the editable sagelib build, which
# re-runs ninja on import, compiles against its numpy headers.
FROM scratch
COPY --from=build / /
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    LC_ALL=C.UTF-8 LANG=C.UTF-8 SHELL=/bin/bash HOME=/home/sage
USER sage
WORKDIR /home/sage/sage
CMD ["/bin/bash"]
