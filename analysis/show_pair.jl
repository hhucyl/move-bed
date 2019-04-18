using HDF5
using PyPlot


PyPlot.clf()
prefix = "/home/pzhang/chen/move-bed/"
middle = "test_mvbed_"
i = 1
name = string(prefix,middle,lpad(i,4,"0"),".h5")
println(name)
f = h5open(name,"r")
Np = 100
pos = read(f["Pposition"])
Pos = [vcat(pos[3*(i-1)+1],pos[3*(i-1)+2],pos[3*(i-1)+3]) for i=1:Np*2]
Ptag = read(f["PTag"])
# flag = true
# IJulia.clear_output(flag)
for ii = 1:2*Np
    if Ptag[ii]>0
        display(scatter!((Pos[ii,1][1],Pos[ii,1][2]),marker=(5,0.2)))

    end
end
