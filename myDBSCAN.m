function [cluster, noise, clusterID] = myDBSCAN(P,eeps,nPts)

clusterID = 0;
n = size(P,1); %nr. total de puncte analizate
cluster = zeros(n,1); %pentru stocarea indexului clusterului corespunzator fiecarui punct
isVisited = zeros(n,1); %pentru a tine cont de punctele analizate
noise = zeros(n,1);

%se parcurg toate punctele ce nu au fost anterior analizate si se verifica numarul de vecini pe o raza eeps
for i=1:n
    if isVisited(i)==0
        isVisited(i) = 1;
        dist = pdist2(P(i,:),P);
        vecini = find(dist<=eeps);
        if numel(vecini)<nPts
            noise(i)=1; %zgomot - daca nr. de vecini este mai mic decat pragul stabilit nPts
        else
            clusterID = clusterID+1; %extind clusterul
            expandCluster(i,vecini,clusterID);
        end
    end
end

function expandCluster(i,vecini,clusterID)
    cluster(i) = clusterID;
    k=1;
    while true
        vecin_curent = vecini(k);
        if isVisited(vecin_curent)==0
            isVisited(vecin_curent) = 1;
            dist_vecini = pdist2(P(vecin_curent,:),P);
            vecini_noi = find(dist_vecini<=eeps);
            if numel(vecini_noi)>=nPts
                vecini = [vecini, vecini_noi];
            end
        end
        if cluster(vecin_curent)==0
            cluster(vecin_curent) = clusterID;
        end
        
        k=k+1;
        if k>numel(vecini)
            break;
        end
    end
end

end





