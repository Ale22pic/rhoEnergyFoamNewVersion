       e = rhoE/rho - 0.5*magSqr(U);
       const labelUList& own   = mesh.owner();
       const labelUList& neigh = mesh.neighbour();
//
//     Loop on all cells
       forAll(neigh,iface)
       {
         //  Left and Right indeces of cell iface
	 label me = own[iface];
	 label my_nei  = neigh[iface];

         scalar hl = H[me] ;
         scalar hr = H[my_nei] ;
         scalar pl = p[me] ;
         scalar pr = p[my_nei] ;
         scalar ul  = U[me]    & (mesh.Sf()[iface]/(mesh.magSf()[iface])) ;
         scalar ur  = U[my_nei]& (mesh.Sf()[iface]/(mesh.magSf()[iface])) ;
	 scalar uul = U[me].component(0)  ;
         scalar uur = U[my_nei].component(0);
         scalar vvl = U[me].component(1)  ;
         scalar vvr = U[my_nei].component(1);
         scalar wwl = U[me].component(2)  ;
         scalar wwr = U[my_nei].component(2);
         scalar dl = rho[me] ;
         scalar dr = rho[my_nei] ;
         scalar el = e[me] ;
         scalar er = e[my_nei] ;


         scalar cj = 0.25*(dl+dr)*(ul+ur) ;
         scalar ij = 0.5*cj*(el+er) ;
	 scalar kj = cj*(uul*uur + vvl*vvr + wwl*wwr)*0.5 ;
         scalar pj = 0.5*(ul*pr + pl*ur) ;

         enFlux[iface] = (ij + kj + pj)*mesh.magSf()[iface] ;
//       enFlux[iface] = 0.125*(dl+dr)*(hl+hr)*(ul+ur) ; 
       }

