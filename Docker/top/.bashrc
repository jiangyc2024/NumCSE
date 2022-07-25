# convenience aliases to improve working with cmake

alias cdbin='cd /numcse/bin$(numcse_path)'
alias cdbuild='cd /build$(numcse_path)'
alias cdsrc='cd /numcse$(numcse_path)'
alias make='make -j$(nproc)'
alias tidy='cmake -DTIDY=1 /build'
alias untidy='cmake -DTIDY=0 /build'