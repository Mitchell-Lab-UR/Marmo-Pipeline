function tim = find_strobe_time(tagcode,taglet,strobes,tstrobes)
%function tim = find_strobe_time(tagcode,sixlet,strobes)
%   find the time in a ephys continuous stream to match a strobe
% input:
%     tagcode - unique code identifying start (63) or end (62) of trial
%     taglet - sixlet of times (year-2000,mo,day,hour,minute,sec) to find
%     strobes - list of integer strobes from 0 to 63
%     tstrobes - list of ephys timings of the strobes
% output:
%     tim - returns the matching time of strobe, or NaN
   zz = find(strobes == tagcode);
   zN = size(zz,1);
   N = size(strobes,1);
   ctaglet = char( taglet );
   %**********
   match = [];
   for k = 1:zN
      kk = zz(k);
      kN = kk+6;
      if (kN <= N)
         mcode = char( strobes((kk+1):(kk+6))' );
         if (strcmp(ctaglet,mcode))
            match = [match kk]; 
         end
      end  
   end
   %***** check for errors in synching codes and report
   if isempty(match)
       % disp(sprintf('No match found for code %d : ',tagcode));
       % taglet
       tim = NaN;
       return;
   end
   if (size(match,2) > 1)
       % disp(sprintf('Multiple matchs found for code %d : ',tagcode));
       % taglet
       tim = NaN;
       return;       
   end
   %*****************
   tim = tstrobes(match);
   %*************
end

