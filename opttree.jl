using DataStructures

function opttree(id)
    ysol, csol = nothing, nothing

    id == 237 && ( (ysol, csol) = (SortedDict(1=>'D', 2=>'*', 3=>'*', 4=>1, 5=>2, 6=>'P', 7=>'C'), SortedDict(7=>2.0)) )

    ysol, csol
end

# ysol, csol = opttree(237)
# println(ysol)
# println(csol)