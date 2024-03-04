	// Collide particles within each cell (here's the physics!):
	function collide() {
		var viscosity = Number(viscSlider.value);	// kinematic viscosity coefficient in natural units
		var omega = 1 / (3*viscosity + 0.5);		// reciprocal of relaxation time
		for (var y=1; y<ydim-1; y++) {
			for (var x=1; x<xdim-1; x++) {
				var i = x + y*xdim;		// array index for this lattice site
				var thisrho = n0[i] + nN[i] + nS[i] + nE[i] + nW[i] + nNW[i] + nNE[i] + nSW[i] + nSE[i];
				rho[i] = thisrho;
				var thisux = (nE[i] + nNE[i] + nSE[i] - nW[i] - nNW[i] - nSW[i]) / thisrho;
				ux[i] = thisux;
				var thisuy = (nN[i] + nNE[i] + nNW[i] - nS[i] - nSE[i] - nSW[i]) / thisrho;
				uy[i] = thisuy
				var one9thrho = one9th * thisrho;		// pre-compute a bunch of stuff for optimization
				var one36thrho = one36th * thisrho;
				var ux3 = 3 * thisux;
				var uy3 = 3 * thisuy;
				var ux2 = thisux * thisux;
				var uy2 = thisuy * thisuy;
				var uxuy2 = 2 * thisux * thisuy;
				var u2 = ux2 + uy2;
				var u215 = 1.5 * u2;
				n0[i]  += omega * (four9ths*thisrho * (1                        - u215) - n0[i]);
				nE[i]  += omega * (   one9thrho * (1 + ux3       + 4.5*ux2        - u215) - nE[i]);
				nW[i]  += omega * (   one9thrho * (1 - ux3       + 4.5*ux2        - u215) - nW[i]);
				nN[i]  += omega * (   one9thrho * (1 + uy3       + 4.5*uy2        - u215) - nN[i]);
				nS[i]  += omega * (   one9thrho * (1 - uy3       + 4.5*uy2        - u215) - nS[i]);
				nNE[i] += omega * (  one36thrho * (1 + ux3 + uy3 + 4.5*(u2+uxuy2) - u215) - nNE[i]);
				nSE[i] += omega * (  one36thrho * (1 + ux3 - uy3 + 4.5*(u2-uxuy2) - u215) - nSE[i]);
				nNW[i] += omega * (  one36thrho * (1 - ux3 + uy3 + 4.5*(u2-uxuy2) - u215) - nNW[i]);
				nSW[i] += omega * (  one36thrho * (1 - ux3 - uy3 + 4.5*(u2+uxuy2) - u215) - nSW[i]);
			}
		}
		for (var y=1; y<ydim-2; y++) {
			nW[xdim-1+y*xdim] = nW[xdim-2+y*xdim];		// at right end, copy left-flowing densities from next row to the left
			nNW[xdim-1+y*xdim] = nNW[xdim-2+y*xdim];
			nSW[xdim-1+y*xdim] = nSW[xdim-2+y*xdim];
		}
	}

	// Move particles along their directions of motion:
	function stream() {
		barrierCount = 0; barrierxSum = 0; barrierySum = 0;
		barrierFx = 0.0; barrierFy = 0.0;
		for (var y=ydim-2; y>0; y--) {			// first start in NW corner...
			for (var x=1; x<xdim-1; x++) {
				nN[x+y*xdim] = nN[x+(y-1)*xdim];			// move the north-moving particles
				nNW[x+y*xdim] = nNW[x+1+(y-1)*xdim];		// and the northwest-moving particles
			}
		}
		for (var y=ydim-2; y>0; y--) {			// now start in NE corner...
			for (var x=xdim-2; x>0; x--) {
				nE[x+y*xdim] = nE[x-1+y*xdim];			// move the east-moving particles
				nNE[x+y*xdim] = nNE[x-1+(y-1)*xdim];		// and the northeast-moving particles
			}
		}
		for (var y=1; y<ydim-1; y++) {			// now start in SE corner...
			for (var x=xdim-2; x>0; x--) {
				nS[x+y*xdim] = nS[x+(y+1)*xdim];			// move the south-moving particles
				nSE[x+y*xdim] = nSE[x-1+(y+1)*xdim];		// and the southeast-moving particles
			}
		}
		for (var y=1; y<ydim-1; y++) {				// now start in the SW corner...
			for (var x=1; x<xdim-1; x++) {
				nW[x+y*xdim] = nW[x+1+y*xdim];			// move the west-moving particles
				nSW[x+y*xdim] = nSW[x+1+(y+1)*xdim];		// and the southwest-moving particles
			}
		}
		for (var y=1; y<ydim-1; y++) {				// Now handle bounce-back from barriers
			for (var x=1; x<xdim-1; x++) {
				if (barrier[x+y*xdim]) {
					var index = x + y*xdim;
					nE[x+1+y*xdim] = nW[index];
					nW[x-1+y*xdim] = nE[index];
					nN[x+(y+1)*xdim] = nS[index];
					nS[x+(y-1)*xdim] = nN[index];
					nNE[x+1+(y+1)*xdim] = nSW[index];
					nNW[x-1+(y+1)*xdim] = nSE[index];
					nSE[x+1+(y-1)*xdim] = nNW[index];
					nSW[x-1+(y-1)*xdim] = nNE[index];
					// Keep track of stuff needed to plot force vector:
					barrierCount++;
					barrierxSum += x;
					barrierySum += y;
					barrierFx += nE[index] + nNE[index] + nSE[index] - nW[index] - nNW[index] - nSW[index];
					barrierFy += nN[index] + nNE[index] + nNW[index] - nS[index] - nSE[index] - nSW[index];
				}
			}
		}
	}
