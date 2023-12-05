function [decoded_codeword, error_vec] = decode_llr(obj, input_llr_vec, max_iter, min_sum)
    eta = zeros(obj.M, obj.N);
    lasteta = zeros(obj.M, obj.N);
    updated_llr_vec = input_llr_vec;
    error_vec = zeros(max_iter, 1);
    for iter = 1 : max_iter
    
        for i_m = 1 : obj.M
            for i_n1 = 1 : obj.row_weight_vec(i_m) 
                n1 = obj.row_mat(i_m, i_n1);
                if min_sum
                    pr = 100;
                else
                    pr = 1;
                end
                for  i_n2 = 1 : obj.row_weight_vec(i_m)
                    if i_n1 == i_n2
                        continue;
                    end
                    n2 = obj.row_mat(i_m, i_n2);
                    l1 = (updated_llr_vec(n2) - lasteta(i_m, n2));
                    l1 = min(l1, 20);
                    l1 = max(l1, -20);
                    if min_sum
                        pr = sign(pr) * sign(l1) * min(abs(l1), abs(pr));
                    else
                        pr = pr * tanh(l1/2);
                    end
                end
                if min_sum
                     eta(i_m, n1) = pr; 
                else
                    eta(i_m, n1) = 2 * atanh(pr);
                end
    
            end
         end
         
         lasteta = eta;
         
         for i_n = 1 : obj.N
             updated_llr_vec(i_n) = input_llr_vec(i_n);
             for i_m = 1 : obj.col_weight_vec(i_n)
                 m = obj.col_mat(i_n, i_m);
                 updated_llr_vec(i_n) = updated_llr_vec(i_n) + eta(m,i_n);
             end
         end
         
         decoded_codeword = (updated_llr_vec < 0);
         if obj.check_codeword(decoded_codeword)
              return;
         else
             error_vec(iter) = 1;
         end
    end

end


function [b] = check_codeword(obj, x)
   b = 1;
   for i_check = 1 : obj.M
       c = 0;
       for i_n = 1  : obj.row_weight_vec(i_check)
           c = c + x(obj.row_mat(i_check, i_n));
       end
       if mod(c, 2) == 1
           b = 0;
           break;
       end                
   end
end