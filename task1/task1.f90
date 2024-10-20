
! cos(x/10) * arctg(x)
! c = 0.5 , d = 1.5
! N = 20 , n0 = 2 , m = 3
module compute_module
    implicit none
contains

    function compute_func(c,i,h)
        real, intent(in) :: c, h
        integer, intent(in) :: i
        real :: compute_func
        compute_func = cos((c + i*h)/10) * atan(c+i*h)
    end function compute_func

    function compute_scalar(c,i,h,vector,size) 
        ! m - number of coefficients (dimension of the polynomial)
        real, intent(in) :: c, h
        integer, intent(in) :: i, size
        integer :: order_number
        real :: vector(size)
        real :: compute_scalar

        ! calculate the sum a0*x^n + a1*x^(n-1) +... + an
        compute_scalar = 0.0
        do order_number = 1,size
            compute_scalar = compute_scalar + vector(order_number) * ((c + i*h)**(size-order_number))
        end do

    end function compute_scalar
end module compute_module



program call_sub
    use compute_module
    implicit none
    real :: c, d, h, eps , func, F_k , derevative_value ,  lambda_k, delta_k_1, beta_k, module_grad_k, module_grad_k_1 , y_i, sigma_n , derev_k_1,  derev_k, term_1, term_2, term_3
    integer :: N_big , n0, m_number, i, j, kol_iterations , m
    real , allocatable:: v_k(:) , v_k_1(:), x_k_1(:), x_k(:), grad_k(:),  grad_k_1(:)

    c = 0.5
    d = 1.5

    N_big = 20 ! number of split points

    n0 = 2! initial degree of polynomial
    m_number = 3 ! need to construct 3 polynomials with degrees 2,3,4

    eps = 0.0001 ! accuracy of approximation 
 
    h = (d - c)/N_big ! grid step size

    print *
    print *,'           Table of values:'
    do i = 0, N_big
        func = compute_func(c,i,h)
        print "(A, F7.5, A, F7.5)", 'Value cos(x/10)*arctg(x) = ',func,' with x = ', c+i*h
    end do  
    print *

    do m = n0+1, m_number+n0
        ! allocate memory for vectors
        allocate(v_k(m))
        allocate(v_k_1(m))
        allocate(x_k_1(m))
        allocate(x_k(m))
        allocate(grad_k(m))
        allocate(grad_k_1(m))

        if (m==m_number) then
            ! x_0 = 0 , if our first polynomial is n=n0
            do i = 1,m
                x_k(i) = 0.0
            end do
        else
            ! if not the first polynomial then a_0 = 0, a_1 = x_1, a_2 = x_2, ...
            x_k(1) = 0.0
            do i = 2, m
                x_k(i) = x_k_1(i)
            end do
        end if

        delta_k_1 = 1000
        kol_iterations = 1
        
        ! find the gradient  v_0 = grad(F(x(0)))
        do j = 1, m
            derevative_value = 0.0
            do i = 0,N_big
                derevative_value = derevative_value +  (compute_func(c,i,h) - compute_scalar(c,i,h,x_k,m))*((c + i*h)**(m-j))
            end do
            v_k(j) = -2.0 * derevative_value
        end do


        ! find |grad(F(x_0))|
        module_grad_k = 0.0
        do i = 1,m
            module_grad_k = module_grad_k + (v_k(i)**2)
        end do 
        module_grad_k = sqrt(module_grad_k)

        ! v_0 = -grad(F(x(0)))
        do i = 1,m
            v_k(i) = -v_k(i)
        end do

        ! find lambda_k, to find later x_k_1
        term_1 = 0.0
        term_2 = 0.0
        term_3 = 0.0    
        lambda_k = 0.0
        do i = 0,N_big
            term_1 = term_1 + compute_func(c,i,h)*compute_scalar(c,i,h,v_k,m)
            term_2 = term_2 + compute_scalar(c,i,h,v_k,m)*compute_scalar(c,i,h,x_k,m)
            term_3 = term_3 + compute_scalar(c,i,h,v_k,m)*compute_scalar(c,i,h,v_k,m)
        end do
        lambda_k = (term_1 - term_2)/term_3
        

        ! find x_k_1 == x_1
        do i = 1,m
            x_k_1(i) = x_k(i) + lambda_k*v_k(i)
        end do
        ! evaluate the accuracy of approximation to the minimum 
        delta_k_1 = 0.0
        do i=1,m
            delta_k_1 = delta_k_1 + ((x_k_1(i) - x_k(i))**2)
        end do
        delta_k_1 = sqrt(delta_k_1)

        print "(A, I2, A, F8.5, A, F8.5)" , " iteration = ", kol_iterations,  " ; |grad(x_k)| = ", module_grad_k ,  " ; delta_k_1 = ",delta_k_1 
    
        kol_iterations = kol_iterations + 1

        do
            if (delta_k_1 < eps ) exit ! find the answer - x_k_1

            ! reset gradients  
            do i = 1,m
                grad_k(i) = 0.0
                grad_k_1(i) = 0.0
            end do

            ! find grad(F(x_k_1)) and  grad(F(x_k))  
            do j = 1, m
                derev_k_1 = 0.0
                derev_k = 0.0
                do i = 0, N_big
                    derev_k_1 = derev_k_1 +  (compute_func(c,i,h) - compute_scalar(c,i,h,x_k_1,m)) * ((c + i*h)**(m-j))
                    derev_k = derev_k +  (compute_func(c,i,h) - compute_scalar(c,i,h,x_k,m)) * ((c + i*h)**(m-j))
                end do
                grad_k_1(j) = -2.0 * derev_k_1
                grad_k(j) = -2.0 * derev_k
            end do


            ! find |grad(F(x_k_1))|^2
            module_grad_k_1 = 0.0
            do i = 1,m
                module_grad_k_1 = module_grad_k_1 + (grad_k_1(i)**2)
            end do 
            module_grad_k_1 = sqrt(module_grad_k_1)
            
            ! find |grad(F(x_k))|^2
            module_grad_k = 0.0
            do i = 1,m
                module_grad_k = module_grad_k + (grad_k(i)**2)
            end do 
            module_grad_k = sqrt(module_grad_k)

            beta_k = (module_grad_k_1**2)/(module_grad_k**2) ! because of  |grad(F(x_k_1))|^2 / |grad(F(x_k))|^2


            ! find v_k_1 
            do i = 1,m
                v_k_1(i) = (-grad_k_1(i)) + (beta_k*v_k(i)) 
            end do 

            ! make substitutions to find the next one lambda_k
            x_k = x_k_1
            v_k = v_k_1

            ! find lambda_k, to find later x_k_1
            term_1 = 0.0
            term_2 = 0.0
            term_3 = 0.0    
            lambda_k = 0.0
            do i = 0,N_big
                term_1 = term_1 + compute_func(c,i,h)*compute_scalar(c,i,h,v_k,m)
                term_2 = term_2 + compute_scalar(c,i,h,v_k,m)*compute_scalar(c,i,h,x_k,m)
                term_3 = term_3 + compute_scalar(c,i,h,v_k,m)*compute_scalar(c,i,h,v_k,m)
            end do
            lambda_k = (term_1 - term_2)/term_3


            ! find x_k_1
            do i = 1,m
                x_k_1(i) = x_k(i) + (lambda_k*v_k(i))
            end do


            ! evaluate the accuracy of approximation to the minimum:
            delta_k_1 = 0.0
            do i = 1, m
                delta_k_1 = delta_k_1 + ( (x_k_1(i) - x_k(i)) **2)
            end do
            delta_k_1 = sqrt(delta_k_1)

            print "(A, I2, A, F8.5, A, F8.5)", " iteration = ", kol_iterations,  " ; |grad(x_k)| = ", module_grad_k , " ; delta_k_1 = ",delta_k_1 
            kol_iterations = kol_iterations + 1

        end do

        print *
        print *,"Answer: x_k_1 = ", x_k_1


        do i = 0,N_big
            y_i = 0.0
            y_i = compute_scalar(c,i,h,x_k,m)
            print "(A, I2, A, F8.5)", " y_",i, " =",  y_i
        end do

        F_k = 0.0 
        do i= 0,N_big
            F_k = F_k + (compute_func(c,i,h) - compute_scalar(c,i,h,x_k,m))**2
        end do

        sigma_n = sqrt(F_k/(N_big+1))
        print "(A, F7.5)"," sigma_n = ", sigma_n
        print *
        print *
        
        ! freeing up memory for vectors
        deallocate(v_k)
        deallocate(v_k_1)
        deallocate(x_k)
        deallocate(x_k_1)
        deallocate(grad_k)
        deallocate(grad_k_1)
    end do
end program call_sub