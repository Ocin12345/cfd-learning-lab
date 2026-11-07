"""Explicit gravitational work for the pressure-based ideal-gas configuration tested here."""
def install_gravity_work(st,g=9.81):
    expression=f'-{g:.17g} [m s^-2] * Density * Velocity.y'
    sources=st.setup.cell_zone_conditions.fluid['fluid'].sources
    sources.enable=True
    sources.terms['energy'].resize(size=1)
    sources.terms['energy'][0].option='value'
    sources.terms['energy'][0].value=expression
    state=sources.get_state()
    assert state['enable'] is True
    assert len(state['terms']['energy'])==1
    assert state['terms']['energy'][0]['value']==expression
    return dict(expression=expression,units='W/m^3',source='rho g dot U for gravity (0,-g)',readback=state)
