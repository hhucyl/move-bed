function kg = Ghost(k,Np,ptag)
    kg = k;
    if(ptag(k(1))>0&&ptag(k(2))<0)
        kg = [k(1)+Np,k(2)];
    else
        kg = [k(1),k(2)+Np];
    end
end