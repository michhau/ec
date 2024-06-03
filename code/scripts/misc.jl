

#################################################################
#plot log wind profile

using LsqFit

h = [0.3, 1, 2, 2.85]#[1, 2, 2.85]
means = [mean(evaldf5.u), mean(evaldf2.u), mean(evaldf3.u), mean(evaldf4.u)]#[mean(evaldf2.u), mean(evaldf3.u), mean(evaldf4.u)]#[mean(evaldf5.u), mean(evaldf2.u), mean(evaldf3.u), mean(evaldf4.u), mean(evaldf6.u)]

logwindprof(z, p) = p[1] / 0.4 * (log.(z / p[2])) #.+ p[3]) #p[1]=u⋆, p[2]=z₀
p0 = [0.2, 1e-3] #initial parameters
logwindfit = curve_fit(logwindprof, h, means, p0)
logwindparam = logwindfit.param
zdata = collect(0.01:0.01:5)

fig = PyPlot.figure()
ax = fig.add_subplot(111)
ax.plot(mean0_3m, 0.3, ".")
ax.plot(mean1m, 1, ".")
ax.plot(mean2m, 2, ".")
ax.plot(mean3m, 3, ".")
ax.plot(mean5m, 5, ".")
ax.plot(logwindprof(zdata, logwindparam), zdata)
ax.plot(log(5/5e-3)^(-1)*mean(evaldf6.u).*log.(zdata./5e-3), zdata)
ax.set_xlabel(L"\overline{u}~\mathrm{[m~s^{-1}]}")
ax.set_ylabel(L"h~\mathrm{[m]}")
ax.grid()
ax.legend()

#########################################################################
#plotting wind direction histograms
periodstart = DateTime(2021,06,08,10,00,00)
periodend = DateTime(2021,06,11,08,17,15)

ec1 = copy(evaldf4)
ec2tmp = copy(t2irgdb)

ec1 = ec1tmp[periodstart .<= ec1tmp.time .<= periodend,:]
ec2 = ec2tmp[periodstart .<= ec2tmp.time .<= periodend,:]
turb.drdf!(ec1)
turb.drdf!(ec2)
#ec3 = tjkdf[periodstart .<= tjkdf.time .<= periodend,:]
#tjkmeteo = tjkmeteodata[]

winddir1 = turb.winddir(ec1)
winddir2 = turb.winddir(ec2)
#winddir3 = turb.winddir(tjkdf)

winddirfig = PyPlot.figure()
dirax = winddirfig.add_subplot(111)
dirax.set_xlabel("wind direction [deg in instrument coord.]")
dirax.set_title("Wind direction histogram 8.6.-11.6. - before DR")
dirax.hist(winddir1.α, bins=collect(0:0.5:360), density=true, label="T2UCSAT", alpha=0.5)
#dirax.hist(winddir2.α, bins=collect(0:0.5:360), density=true, label="T2IRG", alpha=0.5)
#dirax.hist(tjkmeteo.wind_mean_vector_direction, bins=collect(0:2:360), density=true, label="TJK_prop", alpha=0.5)
#dirax.hist(winddir3.α, bins=collect(0:0.5:360), density=true, label="TJK_3DUS", alpha=0.5)
dirax.legend()
PyPlot.show()

PyPlot.plot(winddir1.time[1:10000], abs.(winddir1.α[1:10000]-winddir2.α[1:10000]))
