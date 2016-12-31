function [ output_args ] = Hough( D, Cores )
%HOUGH Shows the Hough space for a given data set D
% close all
%figure
%hold on
StepSize=1024;
    th=[-pi/2:pi/StepSize:pi/2];
    
    n = max(size(D));
    
    R = @(a,b,T) a*cos(T)+b*sin(T);
    MaxR = 0;
    parfor i=1:n
    for i=1:n
        
        X=D(i,1); Y=D(i,2);
        
        RList(i) = abs(max(R(X,Y,th)));
        
    end
    MaxR=max(RList);
    
Zy = -MaxR:2*MaxR/StepSize:MaxR;
Zx = th;
Z=zeros(max(size(Zx)), max(size(Zy)));

spmd
    for i=labindex:Cores:n
        
        if mod(i,1000)==0
            disp(100*i/n)
        end
        
        X=D(i,1); Y=D(i,2);
         
        for j=1:max(size(th))
            %find nearest point in Z to the given and accumulate
            
            t0 = th(j);
            r0 = R(X,Y,t0);

            Diffs = abs(Zx - t0);
            Min = min(Diffs);  Min = find(Diffs==Min);
            Zxloc = Min(1);
            
            Diffs = abs(Zy - r0);
            Min = min(Diffs);  Min = find(Diffs==Min);
            Zyloc = Min(1);
            
            Z(Zxloc,Zyloc) = Z(Zxloc,Zyloc) + 1;
            
        end     
    end
end

Z2 = Z{1};
for i=2:Cores
    Z2 = Z2 + Z{i};
end
Z=Z2;


%     for i=1:n
%         X=D(i,1); Y=D(i,2);
%          
%         for j=1:max(size(th))
%             %find nearest point in Z to the given and accumulate
%             
%             t0 = th(j);
%             r0 = R(X,Y,t0);
% 
%             Diffs = abs(Zx - t0);
%             Min = min(Diffs);  Min = find(Diffs==Min);
%             Zxloc = Min(1);
%             
%             Diffs = abs(Zy - r0);
%             Min = min(Diffs);  Min = find(Diffs==Min);
%             Zyloc = Min(1);
%             
%             Z(Zxloc,Zyloc) = Z(Zxloc,Zyloc) + 1;
%             
%         end
%         
%     end






%     Zx(1)=[]; Zx(end)=[];
%     Zy(1)=[]; Zy(end)=[];
%     Z(1,:)=[];Z(end,:)=[];
%     Z(:,1)=[];Z(:,end)=[];
    imagesc(Zx,Zy,Z');
    
    xlabel('Angle');ylabel('Radius');
    title('Accumulated Hough Space');
%    
%    figure; contour(Zx,Zy,Z');
%    
%    xlabel('Angle');ylabel('Radius')
%    title('Accumulated Hough Space Contour')

  %get best line
  Max = max(max(Z));
  [ROW,COL]=find(Z==Max);
  ROW=ROW(1);COL=COL(1);
  
  m = -1*cot( Zx(ROW) )
  b = Zy(COL) / sin( Zx(ROW) )
  g=@(z) m*z+b;
  
%       
%   figure
%   Xs=D(:,1);
%   plot(D(:,1),D(:,2),'.',Xs,g(Xs),'--')
%   xlabel('x');ylabel('y');
%   title(['Data with most probable line: m=',num2str(m),' b=', num2str(b)])

end

