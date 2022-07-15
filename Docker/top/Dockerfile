from nmcse-dev-base

# ==== prepare repository ==========================================================

# preheat cmake build process

copy . ./numcse
run mkdir /build && cd /build \
      && cmake -DTIDY=0 /numcse \
      && rm -r /numcse

copy Docker/top/numcse_path /bin
copy Docker/top/.bashrc /tmp/patch
run cat /tmp/patch >> ~/.bashrc && rm /tmp/patch

workdir /build