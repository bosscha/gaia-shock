## functions to be used mainly by the build reprocessing scripts

function extra(m::GaiaClustering.meta, optim)
    tstart= now()
    rng = MersenneTwister()
    uuid=uuid4(rng)
    m.uuid= uuid

    println("###########################")
    println("## Starting with $(m.votname)")
    println("## Starting at $tstart")
    println("## Id $uuid")

    df , dfcart , dfcartnorm = get_data(m)

    cycle, flag= cycle_extraction_optim(df, dfcart, m, optim)

    tend= now()
    println("## Ending at $tend")
    println("## number of cycle: $cycle , flag:$flag ")
    println("##########################")
    println("##")

    duration= Dates.value(tend-tstart) / (1000*3600)
    durationstr= @sprintf("%3.3f", duration)
    @printf("## %s \n",specialstr("Duration: $durationstr hours","YELLOW"))
    @printf("## %s \n",specialstr("Votable done: $(m.votname)","YELLOW"))
    println("##\n##")
end