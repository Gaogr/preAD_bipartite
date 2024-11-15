%%%%%%%%%%%%by Lyuan Xu
function matA = bip_proj(matGM,thre)

mat1 = matGM;
mat1(matGM<thre) = NaN;
matA = zeros(size(mat1,2),size(mat1,2));

for ii = 1 : size(mat1,2)
    for jj = ii : size(mat1,2)
        if ii == jj
            matA(ii,jj) = NaN;
        else
            tmp1 = mat1(:,ii);
            tmp2 = mat1(:,jj);
            ttmark = ~isnan(tmp1.*tmp2);
            matA(ii,jj) = sum(tmp1(ttmark));
            matA(jj,ii) = sum(tmp2(ttmark));
        end
    end
end

end