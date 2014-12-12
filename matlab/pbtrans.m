function pbfinal=smdp_makepbtrans(numnodes,blmax,pb,nzmax)
  %{
      %inputs:
      numnodes: number of nodes in the network, excluding sink
      blmax: array of length numnodes, blmax{i} is the number of discrete
              battery levels of node i
      pb: a cell array of size 1xnumnodes; each cell pb{i} is the battery
              level transition probability matrix of node i
      nzmax: maximum nonzero values our sparse matrices can have
              (strictly for implementation issues to avoid out of memory
              problems in matlab)
  
      %---debug values---
      nzmax=10000;
      numnodes=3;
      blmax=[2 2 2];
      pb{1}=eye(3);
      pb{2}=[0.5 0.5 0;0 0.5 0.5; 0.4 0.4 0.2];
      pb{3}=[0.1 0.1 0.8;0.1 0.8 0.1; 0.8 0.1 0.1];
  %}
  
  
  %check that inputs are valid
  %check that each (square) matrix in pb has size to blmax(i)+1
  if length(pb)~=numnodes || length(blmax)~=numnodes
      error('blmax or pb matrix size does not match numnodes');
  end
  for i=1:numnodes
      if size(pb{i},1)~=size(pb{i},2)
          error('each pb trans must be square');
      end
      if size(pb{i},1)~=blmax(i)+1
          error('size of matrix %i does not match blmax',i);
      end
  end
  
  %number of all possible B vectors
  %+1 to include the b0 level; so if blmax(i)=2, possible levels are [0,1,2]
  numbvects=prod(blmax+1);
  
  %Note: Matlab currently does not support sparse 3D matrices (only 2D)
  %make 2 sparse (square) matrices, one to store the final transition
  %probabilities and is updated for each node; another as a temporary working
  %matrix for the node being considered
  pbfinal=spalloc(numbvects,numbvects,nzmax);
  pbtemp=spalloc(numbvects,numbvects,nzmax);
  
  for i=1:numnodes %for each node
      %expand each element in pb{i} by repeating it in a square matrix of
      %size expandcell x expandcell
      expandcell=prod(blmax(i+1:end)+1);
      %also define a temporary matrix for tiling
      pbtile=spalloc(size(pb{i},1)*expandcell,size(pb{i},2)*expandcell,nzmax);
      for j=1:size(pb{i},1)
          for k=1:size(pb{i},2)
              %sanity check
              if pb{i}(j,k)<0
                  error('invalid transition probability value');
              end
              if (pb{i}(j,k)>0)
                  pbtile(1+(j-1)*expandcell:j*expandcell,1+(k-1)* ...
                      expandcell:k*expandcell)=pb{i}(j,k);
              end
          end
      end
      
      %have formed the tile, now repeat the tile as many times to form the
      %numbvects x numbvects matrix
      pbtemp=repmat(pbtile,(blmax(i)+1)^(i-1),(blmax(i)+1)^(i-1));
      
      %if first node, just set pbfinal to pbtemp, otherwise do an 
      %element-wise multiplication of each tiled matrix with pbfinal
      if i==1
          pbfinal=pbtemp;
      else
          pbfinal=pbfinal.*pbtemp;
      end
      %reset the temporary matrices
      pbtemp=[];  
      pbtile=[];
      pbtemp=spalloc(numbvects,numbvects,nzmax);
      pbtile=spalloc(numbvects,numbvects,nzmax);    
  end
  
  %sanity check; each row must total 1
  if(any(abs(sum(pbfinal,2)'-1)>0.0001))
      error('each row should add up to 1; must transition to a new state');
  end
end
