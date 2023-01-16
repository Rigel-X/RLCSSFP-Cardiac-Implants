% 2021/12/24, retrospective undersmapling to simulate higher undersampling
% factor, jie xiang @yale mrrc

function [undersam_kspace] = LcssfpUnderSam(full_kspace)
[oversmpl,Phases,Np,Nx]=size(full_kspace);
for i = 1:oversmpl
    for j = 1:Phases
        for p = 1: Np
            if mod(j,2) == 1
                undersam_kspace(i,j,:,:)=full_kspace(i,j,1:4:end,:);
            else
                undersam_kspace(i,j,:,:)=full_kspace(i,j,3:4:end,:);
            end
        end
    end
end
end

