function Signal_Transformation_Demo
% SIGNAL TRANSFORMATION CHAIN DEMONSTRATOR
% Forward chain : shift -> time scale -> reverse -> amplitude
%                 y(t) = A * x(a*t + b),  a = +/-s,  b = -shift
% Inverse chain : given only x(t) and y(t), estimate (A, a, b) by
%                 minimising the squared difference (least squares).

%% ---------------- Time axis and base signals ----------------
t  = linspace(-10, 10, 1000);          % time vector
dt = t(2) - t(1);                      % sample spacing
x  = [];  y = [];                      % original / final signal (shared)

names = {'Sine','Cosine','Unit Step','Ramp','Gaussian','Impulse'};
sigs  = {@sin, @cos, @(u) double(u >= 0), @(u) max(u,0), ...
         @(u) exp(-u.^2), @(u) double(abs(u) <= dt/2)/dt};

%% ---------------- Figure and control panel ----------------
fig = figure('Name','Signal Transformation Chain Demonstrator', ...
    'NumberTitle','off','Units','pixels','Position',[50 50 1400 780], ...
    'Color',[0.94 0.96 0.98]);

pnl = uipanel(fig,'Units','normalized','Position',[0.01 0.915 0.98 0.075], ...
    'BackgroundColor',[0.86 0.90 0.95],'Title','');

% helper: create a control inside the panel using normalized coordinates
mk = @(style,str,x0,w,y0,h) uicontrol(pnl,'Style',style,'String',str, ...
    'Units','normalized','Position',[x0 y0 w h],'FontSize',10, ...
    'BackgroundColor',get(pnl,'BackgroundColor'));

mk('text','Signal',      0.010,0.050,0.20,0.45);
sigMenu  = mk('popupmenu',names, 0.060,0.100,0.20,0.60);
mk('text','Shift',       0.180,0.045,0.20,0.45);
shiftBox = mk('edit','2',        0.230,0.060,0.20,0.60);
mk('text','Time Scale',  0.310,0.075,0.20,0.45);
scaleBox = mk('edit','2',        0.390,0.060,0.20,0.60);
mk('text','Amplitude',   0.470,0.065,0.20,0.45);
ampBox   = mk('edit','5',        0.540,0.060,0.20,0.60);
revBox   = mk('checkbox','Reverse',0.620,0.090,0.20,0.60);
set(revBox,'Value',1);
btnApply = mk('pushbutton','Apply Chain',   0.730,0.120,0.12,0.76);
btnInv   = mk('pushbutton','Inverse Chain', 0.860,0.120,0.12,0.76);
set(btnApply,'Callback',@applyChain,'FontWeight','bold');
set(btnInv,  'Callback',@inverseChain,'FontWeight','bold');

% status line under the panel
status = uicontrol(fig,'Style','text','Units','normalized', ...
    'Position',[0.01 0.870 0.98 0.035],'HorizontalAlignment','left', ...
    'FontSize',10,'BackgroundColor',get(fig,'Color'), ...
    'String','Select a signal and press "Apply Chain".');

%% ---------------- 2 x 3 grid of axes (normalized layout) ----------------
colX = [0.055 0.385 0.715];            % left edge of each column
rowY = [0.510 0.085];                  % bottom edge of top / bottom row
axW  = 0.255;  axH = 0.30;             % axes width / height

ax = gobjects(1,6);
for k = 1:6
    r = ceil(k/3);  c = mod(k-1,3) + 1;
    ax(k) = axes('Parent',fig,'Units','normalized', ...
        'Position',[colX(c) rowY(r) axW axH],'Box','on');
end

titles = {'1. Original','2. Shift','3. Time Scaling', ...
          '4. Reversal','5. Amplitude Scaling','6. Inverse Reconstruction'};
colors = {[0 .45 .74],[.85 .33 .10],[.49 .18 .56], ...
          [0 .5 0],[0 0 0]};           % one colour per stage

applyChain();                          % draw initial result

%% ================= Forward chain =================
    function applyChain(~,~)
        % ---- read and validate inputs ----
        tau = str2double(get(shiftBox,'String'));
        s   = str2double(get(scaleBox,'String'));
        A   = str2double(get(ampBox,'String'));
        a   = 1 - 2*get(revBox,'Value');       % +1 normal, -1 reversed
        if any(~isfinite([tau s A])) || s <= 0
            errordlg('Enter valid numbers (Time Scale must be > 0).');
            return
        end

        % ---- build the original and the intermediate signals ----
        x = sigs{get(sigMenu,'Value')}(t);
        I = @(u) interp1(t,x,u,'linear',0);    % evaluate x(u), 0 outside range
        Z = cell(1,5);
        Z{1} = x;                              % original      x(t)
        Z{2} = I(t - tau);                     % shift         x(t - tau)
        Z{3} = I(s*t - tau);                   % time scaling  x(s*t - tau)
        Z{4} = I(a*s*t - tau);                 % reversal      x(-s*t - tau)
        Z{5} = A * Z{4};                       % amplitude     A*x(...)
        y    = Z{5};                           % final signal

        % ---- plot every stage ----
        for k = 1:5
            cla(ax(k));
            plot(ax(k),t,Z{k},'Color',colors{k},'LineWidth',1.6);
            style(ax(k),titles{k});
        end

        % ---- clear the inverse-reconstruction plot ----
        cla(ax(6)); legend(ax(6),'off');
        style(ax(6),titles{6});

        set(status,'String',sprintf( ...
            'Applied: shift = %g | scale = %g | reverse = %d | amplitude = %g', ...
            tau,s,a<0,A));
    end

%% ================= Inverse chain =================
    function inverseChain(~,~)
        if isempty(y), return, end

        G  = @(q) interp1(t,x,q(1)*t + q(2),'linear',0);  % x(a*t + b)
        Af = @(z) (z*y') / (z*z' + eps);                  % best amplitude (closed form)
        J  = @(q) sum((y - Af(G(q))*G(q)).^2);            % squared error to minimise

        % multi-start search avoids local minima
        opts = optimset('Display','off','MaxIter',2000,'TolX',1e-8,'TolFun',1e-10);
        best = Inf;  qb = [1 0];
        for a0 = [-3 -2 -1 -0.5 0.5 1 2 3]
            for b0 = -4:2:4
                [q,f] = fminsearch(J,[a0 b0],opts);
                if f < best, best = f; qb = q; end
            end
        end

        z    = G(qb);                   % x(a*t + b) with best parameters
        Ahat = Af(z);                   % best amplitude
        yhat = Ahat * z;                % estimated final signal
        rmse = sqrt(mean((y - yhat).^2));

        % ---- plot: original (grey), final (blue), estimate (red dashed) ----
        cla(ax(6)); hold(ax(6),'on');
        plot(ax(6),t,x,   'Color',[.6 .6 .6],'LineWidth',1.4);
        plot(ax(6),t,y,   'b',  'LineWidth',1.8);
        plot(ax(6),t,yhat,'r--','LineWidth',1.8);
        hold(ax(6),'off');
        style(ax(6),titles{6});
        legend(ax(6),'Original Signal','Final Signal','Estimated Signal', ...
            'Location','northeast');

        % ---- report: shift = -b, scale = |a|, reverse if a < 0 ----
        set(status,'String',sprintf(['Found: shift = %.3f | scale = %.3f | ' ...
            'reverse = %d | amplitude = %.3f  ->  y(t) = %.3f*x(%.3f*t %+.3f)' ...
            '  |  RMSE = %.2g'], ...
            -qb(2),abs(qb(1)),qb(1)<0,Ahat,Ahat,qb(1),qb(2),rmse));

        fprintf('\nINVERSE CHAIN RESULT\n');
        fprintf('y(t) = %.4f * x(%.4f*t %+.4f),  RMSE = %.3g\n', ...
            Ahat,qb(1),qb(2),rmse);
    end

%% ================= Common axes styling =================
    function style(h,ttl)
        title(h,ttl,'FontWeight','bold');
        xlabel(h,'Time (t)'); ylabel(h,'Amplitude');
        grid(h,'on'); xlim(h,[-10 10]);
    end
end