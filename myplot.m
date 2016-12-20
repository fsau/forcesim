% Copyright (c) 2016, Franco Sauvisky
% All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
% 3. The name of the author may not be used to endorse or promote products
%    derived from this software without specific prior written permission.

% THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
% IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
% NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
% THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

disp("Calculating..."); fflush(stdout);

fflag = system("gcc sim.c -O3 -lm -o sim && ./sim > orbit1.dat");
if fflag != 0
	return;
end

disp("Done. Processing..."); fflush(stdout);

fps = 30;

fh = fopen("./orbit1.dat");
raw = fscanf(fh, "%f");
fclose(fh);

n = raw(1); % Number of particles
raw = raw(2:end);

data = reshape(raw, [2*n+1, length(raw)/(2*n+1)])';

tmax = max(data(:,1));
normaldata = [(0:1/fps:tmax)', spline(data(:,1), data(:,2:end),0:1/fps:tmax)'];
datalength = size(normaldata, 1);

baseaxis = [min(min(data(:,[2:2:end]))),max(max(data(:,[2:2:end]))),...
min(min(data(:,[3:2:end]))),max(max(data(:,[3:2:end])))];
mycmap = colormap(lines);

cgif = 0;
while cgif != 3
	cgif = menu("Done.", "Plot real-time", "Create GIF file", "Exit");

	if cgif == 1
		if isempty(findall(0,'Type','Figure'))
			figure()
		end

		clf;
		disp("Click on figure to plot..."); fflush(stdout);
		waitforbuttonpress();

		i = 1;
		tic;
		while i < datalength
			scatter(normaldata(i,[2:2:end]),normaldata(i,[3:2:end]),5,mycmap(1:n,:),"filled");
			grid on
			axis(baseaxis,"equal");
		title(["t = " num2str(normaldata(i,1),"%2.2f")]);
		pause(0.001);
		i = round(fps*toc);
		end

		hold on
		for i = 1:n
			plot(data(:,2*i),data(:,2*i+1),"color", mycmap(i,:), "linewidth", 1);
		end
		hold off
	end

	if cgif == 2
		disp("Printing frames..."); fflush(stdout);
		set(0, 'defaultfigurevisible', 'off');
		mkdir output;
		for t = 1:size(normaldata,1)
			scatter(normaldata(t,[2:2:end]),normaldata(t,[3:2:end]),5,mycmap(1:n,:),"filled");
			grid on
			axis(baseaxis,"equal");
			title(["t = " num2str(normaldata(t,1),"%2.2f")]);
			filename=sprintf('./output/frame%d.png',t);
			print("-r0", filename);
		end
		set(0, 'defaultfigurevisible', 'on');
		disp("Done. Converting and creating gif..."); fflush(stdout);
		system("for i in $(ls ./output/*.png); do convert $i $i.gif; rm $i; done");
		system(["gifsicle --loop --delay=" num2str(round(100/fps)) " $(ls -v1 ./output/*.gif) > output.gif"]);
		system("rm -r ./output");
	end
end

disp("Done!")
