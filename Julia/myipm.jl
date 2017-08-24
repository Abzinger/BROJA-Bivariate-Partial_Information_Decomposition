# myipm.jl
module MyIPM
export my_ipm

type MyIPM_Parameters
    t_0 :: Float64
    μ   :: Float64
end

type MyIPM_Results{FLOAT}
end

function my_ipm{FLOAT}(e::My_Eval, param::MyIPM_Parameters, dummy::FLOAT) :: MyIPM_Results{FLOAT}
    oldprc = precision(BigFloat)
    setprecision(BigFloat,e.bigfloat_nbits)

    r = ipm(e, param, dummy)

    setprecision(BigFloat,oldprc)

    return res;
end



function ipm{FLOAT}(e::My_Eval, param::MyIPM_Parameters, dummy::FLOAT) :: MyIPM_Results{FLOAT}
    local prjdata :: Projector_DATA{FLOAT}
    setup_prj!(prjdata, e,param)

    local q   = Vector{FLOAT}( Ones(e.n)  )  # current feasible solution
    local ∇   = Vector{FLOAT}( zeros(e.n) )  # gradient
    local pr∇ = Vector{FLOAT}( zeros(e.n) )  # projected gradient
    local x   = Vector{FLOAT}( zeros(e.n) )  # q + pr∇
    local δ   = FLOAT(1)/(e.n_x*100)         # distance from boundary
    local η   = δ                            # step length



    local terminate = false
    local t :: FLOAT = param.t_0

    while ! terminate # outer loop

        inner_keep_going = true
        while inner_keep_going # inner loop
            
        end # while inner loop

        # Check termination
        if ....
            terminate = true
        else
            t *= param.μ
        end
    end # outer loop



    return results;
end
