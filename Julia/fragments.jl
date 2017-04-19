type Set_Data{T1,T2,T3}
    X :: Vector{T1}
    Y :: Vector{T2}
    Z :: Vector{T3}
    var_idx :: Dict{Tuple{Int32,Int32,Int32},Int32}
    Set_Data() = new([],[],[],Dict{Tuple{T1,T2,T3},Int32}())
end

to_yz{I<:Integer              }(d::Set_Data, y::I, z::I)         :: Int32    = Int32( y*d.length(Z) + z )
varidx{I₁<:Integer,I₂<:Integer}(d::Set_Data, x::I₁,yz::I₂)       :: Int32    = get( d.var_idx, (x,yz), 0)

function Set_Data{T1,T2,T3,FLOAT}(pdf::Dict{Tuple{T1,T2,T3},FLOAT})
    local var_list::Vector{ Pair{Tuple{Int32,Int32},Int32} }
    local X_set::Set{T1}
    local Y_set::Set{T2}
    local Z_set::Set{T3}

    for (xyz,val) ∈ pdf
        if val > 0
            x,y,z = xyz
            push!(X_set,x)
            push!(Y_set,y)
            push!(Z_set,z)
        end
    end

    d = Set_Data{T1,T2,T3}
    d.X = Vector{T1}( collect(X_set) )
    d.Y = Vector{T1}( collect(Y_set) )
    d.Z = Vector{T1}( collect(Z_set) )

    CONTINUE HERE
    m_XY::Array{FLOAT} = [  sum(view(pdf[x,y,:]))   for x ∈ d.X, y ∈ d.Y  ]
    m_XZ::Array{FLOAT} = [  sum(view(pdf[x,:,z]))   for x ∈ d.X, z ∈ d.Z  ]


    ...................................................................................................


    local counter::Int = 0
            counter += 1




    sizehint!(d.var_idx, counter)

    sort!(d.X)
    sort!(d.Y)
    sort!(d.Z)

    # We want an X x (YxZ) in column major:
    for i_y,y in enumerate(d.Y):
        for i_z,z in enumerate(d.Z):
            for i_x,x in enumerate(d.X):
                if get(pdf,(x,y,z),-1) > 0
                    push!(var_list, (Int32(i_x),to_yz(i_y,i_z)) => Int32(1+length(var_list)) )
                end
            end
        end
    end

    d.var_idx = Dict{Tuple{Int32,Int32},Int32}( var_list )

    return d
    ;
end #^ Set_Data() [constructor]
