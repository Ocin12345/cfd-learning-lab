"""Read coupled linear histories without conflating initial and corrected residuals."""
import math

def parse_blocks(text,start,end):
 out=[];active=None;tol=None
 for line in text.splitlines():
  if line.strip()=='pressure coupled equations:':
   if active is not None:raise ValueError('Unclosed coupled block')
   active=[];tol=None;continue
  parts=line.split()
  if active is None:continue
  if parts and parts[0]=='tol.':
   try:tol=[float(x) for x in parts[1:]]
   except ValueError as e:raise ValueError('Invalid tolerance row') from e
   if len(tol)!=4 or not all(math.isfinite(x) and x>=0 for x in tol):raise ValueError('Invalid tolerance row')
   continue
  if not parts or not parts[0].isdigit():continue
  if len(parts)==5:
   try:vals=[float(x) for x in parts[1:]]
   except ValueError as e:raise ValueError('Invalid linear residual') from e
   if int(parts[0])!=len(active) or not all(math.isfinite(x) and x>=0 for x in vals):raise ValueError('Invalid/missing inner cycle')
   active.append(vals)
  elif len(parts)>=6 and ':' in parts[5]:
   try:vals=[float(x) for x in parts[1:5]]
   except ValueError as e:raise ValueError('Invalid outer residual') from e
   if len(active)<2 or tol is None or not all(math.isfinite(x) and x>=0 for x in vals):raise ValueError('Incomplete/nonfinite coupled block')
   out.append(dict(iteration=int(parts[0]),initial=active[0],terminal=active[-1],cycles=len(active)-1,outer_printed=vals,absolute_tolerance_printed=tol));active=None
 if active is not None:raise ValueError('Unclosed coupled block at EOF')
 if [r['iteration'] for r in out]!=list(range(start,end+1)):raise ValueError('Missing/duplicate/out-of-range outer iteration')
 return out
