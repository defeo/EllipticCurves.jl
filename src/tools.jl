######################################################################
# tools.jl: useful tools for EllipticCurves
######################################################################

function (P::GenMPoly)(x, y, z)
	return evaluate(P, [x, y, z])
end






