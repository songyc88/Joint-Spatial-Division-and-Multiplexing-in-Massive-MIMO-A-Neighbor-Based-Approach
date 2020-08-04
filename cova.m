function R=cova(the,del,N)
D=1/2;
R=zeros(N,N);
for ii=1:N
    for jj=1:N
        if ii>=jj
        f=@(x)exp(-i*2*pi*D*(ii-jj)*sin(x));
        R(ii,jj)=(180/pi)*(1/del)*quad(f,(the-del/2)*pi/180,(the+del/2)*pi/180);
        end
    end
end
for ii=1:N
    for jj=1:N
        if ii<jj
            R(ii,jj)=R(jj,ii)';
        end
    end
end

        


    