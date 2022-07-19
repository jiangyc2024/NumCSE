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

copy Docker/top/entrypoint.sh /entry
entrypoint /entry

workdir /build