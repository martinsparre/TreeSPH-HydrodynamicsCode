; density test at t=0.

pro A
readcol,"rhoh.txt", rho, h
readcol,"rho0_gadget.txt",ID,rhoG,hG,format="int,float,float"

deviations = fltarr(1472)
deviationsh = fltarr(1472)

equation=fltarr(1472)
equationG=fltarr(1472)

for i=0,1471 do begin
	for j=0,1471 do begin
		if (ID[j] eq i+1) then begin
			deviations[i]=(rhoG[j]-rho[i])/rhoG[j]
			deviationsh[i]=(hG[j]-h[i])/h[j]			
		;	print, 4.*3.14159265359/(3.*50*0.000679) * hG[j]^3 * rhoG[j],4.*3.14159265359/(3.*50*0.000679) * h[i]^3 * rho[i]
		endif
	endfor

	equationG[i]=4.*3.14159265359/(3.*50*0.000679) * hG[i]^3 * rhoG[i]
	equation[i]=4.*3.14159265359/(3.*50*0.000679) * h[i]^3 * rho[i]
endfor

print, "Equation test (4.*3.14159265359/(3.*50*0.000679) * h^3 * rho) - gadget:"
print, "mean=",mean(equationG)
print, "stddev=",stdev(equationG)
print, "max=",max(equationG)
print, "min=",min(equationG)
print,""
print, "Equation test - My own code:"
print, "mean=",mean(equation)
print, "stddev=",stdev(equation)
print, "max=",max(equation)
print, "min=",min(equation)

print,""
print, "Density deviations"
print, "mean=",mean(deviations)
print, "stddev=",stdev(deviations)
print, "max=",max(deviations)
print, "min=",min(deviations)

print,""
print, "h deviations"
print, "mean=",mean(deviationsh)
print, "stddev=",stdev(deviationsh)
print, "max=",max(deviationsh)
print, "min=",min(deviationsh)


end
