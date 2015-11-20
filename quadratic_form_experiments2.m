%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% look only at paths in zero force field, but with
%% non-zero initial velocity v0
%%
%% there are two phases: (1) velocity-dominated phase, 
%% where the path is disturbed by the initial velocity.
%% (2) acceleration dominated phase, where the initial
%% velocity has been compensated by acceleration,
%% and where the complete space is reachable
%%
%% we concentrate here only on (1), and try to find its
%% boundary, so that we have points which are phase-(1)
%% reachable, and points which are phase-(2) reachable
%%
%% phase-(1) reachability should be similar to force 
%% influence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ae = 0*[2.2;1.2;0.0];
aenorm = ae/norm(ae);

qzp = rot2q([0;0;1],+pi/2);
qzm = rot2q([0;0;1],-pi/2);

amax = 0.5;
v0=1.1*[0.2;-1.2;0.0];
v0norm = v0/norm(v0);
delta = @(x,y) ((x*ae(1)+y*ae(2))**2)-((x*x+y*y))*(ae'*ae - amax*amax);
omega = @(x,y) (-(x*x+y*y)*(ae'*ae - amax*amax));

dk = -ae'*ae+amax*amax;
L=10;
N=100;
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
B = zeros(N,N);
T = norm(v0)/amax
B(1,1)=1.0;
%for i=1:N
%        for j=1:N
%                for t = 0:0.1:T
%                        p = [X(i);Y(j);0.0]-0.5*ae*t*t-v0*t;
%                        dd = sqrt((p'*p));
%                        if dd <= (0.5*amax*t*t);
%                                B(i,j)=1.0;
%                                break
%                        end
%                end
%        end
%end

imagesc(X,Y,B);
hold on;
plot(0,0,'w','markersize',N/10)
hold on;
plot([0,ae(2)],[0,ae(1)],'w');
hold on;
plot([0,v0(2)],[0,v0(1)],'g');
hold on;

xlabel('y');
ylabel('x');
set(gca,'YDir','normal')

for angle=0.0:pi/1000:2*pi
        ai = amax*ones(size(v0));
        qz = rot2q([0;0;1],angle);
        aq = quaternion(v0norm(1),v0norm(2),v0norm(3));
        [aoff,tmpangle] = q2rot(qz*aq*conj(qz));
        angle
        if aoff*v0 < 0
                T = norm(v0)/norm(aoff*v0);
                if T < 10
                        %t = 0:0.1:T;
                        %xx = v0*t+0.5*amax*aoff'*(t.^2);
                        %plot(xx(2,:),xx(1,:),'w','linewidth',1);
                        hold on;
                        t = T;
                        xx = v0*t+0.5*amax*aoff'*(t.^2);
                        plot(xx(2,:),xx(1,:),'w','linewidth',4);
                        %pause(0.05);
                end
        end
end
k=0;
%angle = pi/2+k*pi/40;
%angle=pi/2;
%angle=3*pi/2-asin(amax/norm(ae));
%qz = rot2q([0;0;1],angle);
%aq = quaternion(aenorm(1),aenorm(2),aenorm(3));
%[aoff,tmpangle] = q2rot(qz*aq*conj(qz));
%xx = 0.5*ae*(t.^2)+v0*t+0.5*amax*aoff'*(t.^2);
%plot(xx(2,:),xx(1,:),'w','linewidth',2);
%hold on;

pause
