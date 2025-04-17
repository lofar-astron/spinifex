from __future__ import annotations

import tempfile
from pathlib import Path

import pytest
from aiohttp import web
from spinifex.download import download_file_http

EXAMPLE_CONTENT = b"Example content for testing."


@pytest.fixture
async def local_file_server(tmp_path):
    tmp_path = Path(tmp_path)
    file_path = tmp_path / "test_file.dat"
    file_path.write_bytes(EXAMPLE_CONTENT)

    path_list: list[Path] = [file_path]
    file_paths = {path.name: path for path in path_list}

    async def generic_handler(request):
        filename = request.match_info["filename"]
        if filename in file_paths:
            return web.FileResponse(path=file_paths[filename])
        return web.Response(status=404, text="File not found")

    def create_handler(filename):
        async def handler(request):
            return await generic_handler(
                request.clone(match_info={"filename": filename})
            )

        return handler

    app = web.Application()
    # Add individual routes for backward compatibility
    for file_path in path_list:
        filename = file_path.name
        app.router.add_get(f"/{filename}", create_handler(filename))

    # Add a catch-all route for any file
    app.router.add_get("/{filename}", generic_handler)

    runner = web.AppRunner(app)
    await runner.setup()
    site = web.TCPSite(runner, "localhost", 0)
    await site.start()

    port = site._server.sockets[0].getsockname()[1]
    url = f"http://localhost:{port}"

    yield url

    await runner.cleanup()


@pytest.mark.allow_hosts(["127.0.0.1", "::1"])
@pytest.mark.asyncio
async def test_download_file_http(local_file_server):
    file_name = "test_file.dat"
    url = f"{local_file_server}/{file_name}"
    with tempfile.TemporaryDirectory() as tmpdir:
        output_file = Path(tmpdir) / "downloaded.dat"
        await download_file_http(url, output_file)

        # Assert file was written with expected content
        assert output_file.exists()
        assert output_file.read_bytes() == EXAMPLE_CONTENT
