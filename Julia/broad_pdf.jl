#broad_pdf.jl
using JSON

function create_broad_pdf(n_x::Int, n_y::Int, n_z::Int, dependence::Float64=0.) :: Dict{Tuple{Int,Int,Int},BigFloat}
    @assert n_x ≥ 2 "create_broad_pdf(): n_x >= 2 needed"
    @assert n_y ≥ 2 "create_broad_pdf(): n_y >= 2 needed"
    @assert n_z ≥ 2 "create_broad_pdf(): n_z >= 2 needed"
    @assert 0. ≤ dependence ≤ 1. "create_broad_pdf(): indep in [0,1] needed"
    pdf = Dict{Tuple{Int,Int,Int},BigFloat}()
    for y = 1:n_y
        for z = 1:n_z
            for x = 1:n_x
                pdf[ (x,y,z) ] = ( 1 + dependence*(rand()-.5) )/( n_x*n_y*n_z )
            end # for z
        end # for z
    end # for y
    return pdf
end # create_broad_pdf()

function broad_pdf_to_file(filename::String, n_x::Int, n_y::Int, n_z::Int, dependence::Float64=0.) :: Void
    open(filename,"w") do file
        pdf = create_broad_pdf(n_x,n_y,n_z,dependence)
        pdflist = []
        for k in keys(pdf)
            (x,y,z) = k
            e = [[x,y,z],pdf[(x,y,z)]]
            push!(pdflist, e )
        end
        JSON.print(file,pdflist)
    end
    return nothing
end # broad_pdf_to_file()

function loop_broad_pdfs_to_files() :: Void
    for n_x = 2:5
        for n_y = n_x : 5*n_x
            n_z = n_y
            filename = "broad_indep--"*string(n_x)*"-"*string(n_y)*"-"*string(n_z)*".dens"
            broad_pdf_to_file(filename,n_x,n_y,n_z)
            filename = "broad_dep--"*string(n_x)*"-"*string(n_y)*"-"*string(n_z)*".dens"
            broad_pdf_to_file(filename,n_x,n_y,n_z,1.)
        end
    end
    return nothing;
end


;
#EOF
