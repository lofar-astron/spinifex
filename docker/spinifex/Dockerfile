ARG BUILD_ENV=no_copy

FROM python:3.11 AS build_no_copy
ADD ../../requirements.txt .
COPY ../.. /work
RUN rm -r /work/dist | true
RUN python -m pip install --user tox
WORKDIR /work
RUN python -m tox -e build

FROM python:3.11 AS build_copy
COPY dist /work/dist

FROM build_${BUILD_ENV} AS build

FROM python:3.11-slim
COPY --from=build /work/dist /dist
RUN python -m pip install /dist/*.whl
