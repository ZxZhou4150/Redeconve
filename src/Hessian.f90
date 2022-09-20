subroutine Hessian_f90(ref,ngenes,ncells,r,hp) bind(C, name="Hessian_f90_")
    INTEGER ngenes,ncells
    REAL(8) ref(ncells,ngenes)
    REAL(8) r(ncells,ncells)
    INTEGER hp

    REAL(8) G(ncells,ncells)
    INTEGER i,j

    do i=1,ncells
        G(i,i) = sum(ref(:,i)**2)+hp*sum(ref(i,:))
    end do
    do i=1,ncells-1
        do j=i+1,ncells
            G(i,j) = 2*(dot_product(ref(:,i),ref(:,j))-hp*r(i,j))
        end do
    end do
    G = G+transpose(G)
end