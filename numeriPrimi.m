%This function has the aim to verify wheter a number is prime or not
function []= numeriPrimi(x)
if x<=3 && x>=0
    disp('The number is prime');
elseif x<0
    disp('The number is negative, we are considering only positive integer number');
else
    rad=round(sqrt(x));
    rest=mod(x,2);
    k=3;
    if rest==0
        disp('The number is not prime, it is divisible almost for 2');
    else
        while k<=rad
            rest2=mod(x,k);
            if rest2==0
                disp('The number is not prime');
                break
            end
            k=k+1;
        end
    if k>rad
        disp('The number is PRIME');
    end
    end
end
end