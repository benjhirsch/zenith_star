import numpy as np

def round_to(number, sig_figs):
    sig_figs -= 1
    try:
        return round(number, -int(np.floor(np.log10(number)) - sig_figs))
    except:
        0

def delta_ra(ra, zenith_ra, zenith_dec):
    #small angle approximation of gnomonic projection to account for high declination
    return (ra - zenith_ra) * np.cos(zenith_dec)

def delta_dec(dec, zenith_dec):
    return dec - zenith_dec

def bv_to_temp(bv):
    #temperature formula from https://arxiv.org/abs/1201.1809
    try:
        temp = 4600*(1/(.92*bv+1.7) + 1/(.92*bv+.62))
    except:
        temp = None
    return temp

def temp_to_rgb(temp):
    #empirical temperature to rgb formula adapted from https://www.celestialprogramming.com/articles/starColors/ColorHelland.js
    try:
        temp = float(temp)
    except:
        return 255, 255, 255
    
    temp = temp * 0.01
    if temp <= 66:
        r = 255
        g = 99.4708025861 * np.log(temp) - 161.1195681661
        if temp <= 19:
            b = 0
        else:
            b = 138.5177312231 * np.log(temp-10) - 305.0447927307
    else:
        r = 329.698727446 * np.power(temp-60, -0.1332047592)
        g = 288.1221695283 * np.power(temp-60, -0.0755148492)
        b = 255
    
    r = int(min([max([r, 0]), 255]))
    g = int(min([max([g, 0]), 255]))
    b = int(min([max([b, 0]), 255]))

    return r, g, b

def zenith_field(display_stars_list, search_x, search_y, zenith_ra, zenith_dec, brightest):
    try:
        import pandas as pd
        import plotly.express as px
    except ImportError:
        print('Plotting requires pandas and plotly. Run "pip install pandas plotly" from the command line.')
        return

    display_stars_df = pd.DataFrame(display_stars_list)
    
    #lower magnitude = brighter star = larger star marker
    size = (22 - display_stars_df['Mag'])*5

    fig = px.scatter(
        display_stars_df,
        x='delta_RA',
        y='delta_Dec',
        custom_data = ['ID', 'Mag', 'Temp', 'RA', 'Dec', 'Alt', 'Az'],
        size=size,
        opacity=1,
        color='rgb',
        color_discrete_map='identity'
    )

    fig.update_traces(
        hovertemplate = (
            "<b>ID: %{customdata[0]}</b><br>"
            "Visual Mag: %{customdata[1]:.2f}<br>"
            "Temp: %{customdata[2]:,d} K<br>"
            "RA: %{customdata[3]:.2f}°<br>"
            "Dec: %{customdata[4]:.2f}°<br>"
            "Alt: %{customdata[5]:.2f}°<br>"
            "Az: %{customdata[6]:.2f}°<br>"
            "<extra></extra>"
        ),
        marker=dict(line=dict(width=1, color='black'))
    )

    fig.update_layout(
        title='Zenith Field (%s %s)' % (len(display_stars_list), 'brightest stars near the zenith' if brightest else 'closest stars to the zenith'),
        xaxis_title='ΔRA (deg)',
        yaxis_title='ΔDec (deg)',
        xaxis=dict(range=[-search_x, search_x], showgrid=False, zeroline=False),
        yaxis=dict(range=[-search_y, search_y], showgrid=False, zeroline=False, scaleanchor='x'),
        plot_bgcolor='black',
        paper_bgcolor='black',
        font=dict(color='white'), showlegend=False
    )

    fig.add_scatter(x=[0], y=[0], mode='markers', hoverinfo='text', text=['Zenith RA: %.2f°, Zenith Dec: %.2f°' % (zenith_ra, zenith_dec)], marker=dict(color='yellow', size=12, symbol='cross'), name='Zenith')

    return fig
