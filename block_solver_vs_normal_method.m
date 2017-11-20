for n = [600,1200,1800,2400,3000,3600,4200,4800]
    block_solver_method = [];
    normal_methor = [];
    for time = 1:10
        %generate D w y and b
        %The reason for me to divide 100 is because randi only generate integer
        %but with divide 100 means it also include non integer
        random_matrix = randi(10^5,n)/100;
        D = abs(diag(diag(random_matrix)));

        w = randi(10^5)/100;
        y = randi(10^5,[n,1])/100;
        b = randi(10^5,[n,1])/100;
        A = D + y*w*transpose(y);
        
        %time of normal method
        normal_start = tic;
        x_1 = A\b;
        normal_end = toc(normal_start);
        
        block_solver_start = tic;
        %step 1
        y_bar = D\y;
        D_bar = 1 + w*transpose(y)*y_bar;
        b_bar = sqrt(w)*transpose(y_bar)*b;
        
        %step 2
        z = D_bar \ b_bar;
        
        %step3
        x_2 = D\(b-sqrt(w)*z*y);
        block_solver_end = toc(block_solver_start);
        %put to list
        normal_method(time) = normal_end;
        block_solver_method(time) = block_solver_end;
    end

    block_solver_time = mean(block_solver_method);
    normal_time = mean(normal_method);
    disp('when n = ')
    disp(n)
    disp('normal method use: ') 
    disp(normal_time)
    disp('block_solver method use: ')
    disp(block_solver_time)
end