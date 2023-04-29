## aliases to run extra packge with docker
alias getgaia='docker run -v $(pwd):/data bosscha/gaia /soft/julia-1.6.7/bin/julia /run/scripts/getgaia.jl'
alias extra='docker run -v $(pwd):/data bosscha/gaia /soft/julia-1.6.7/bin/julia /run/scripts/extra.jl'
alias build='docker run -v $(pwd):/data bosscha/gaia /soft/julia-1.6.7/bin/julia /run/scripts/build.jl'