from numcse_dev_base

# ==== prepare repository ==========================================================

# preheat cmake build process

copy . /numcse
run mkdir /build && cd /build \
       && cmake -DTIDY=0 -DDOCKER=1 /numcse \
       && mkdir -p /bin_cache \
       && cp -r /numcse/bin/* /bin_cache \
       && rm -r /numcse

copy Docker/top/numcse_path /bin
copy Docker/top/.bashrc /tmp/patch
run cat /tmp/patch >> ~/.bashrc && rm /tmp/patch

# TODO: move into base Dockerfile
# install VS Code (code-server)
RUN curl -fsSL https://code-server.dev/install.sh | sh -s -- --version 4.14.1

# install VS Code extensions
RUN code-server \
       --install-extension notskm.clang-tidy \
       --install-extension ms-vscode.cmake-tools \
       --install-extension twxs.cmake \
       --install-extension ms-vscode.cpptools-themes

copy Docker/top/entrypoint.sh /entry
entrypoint /entry

workdir /build