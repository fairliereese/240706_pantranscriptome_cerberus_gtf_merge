import seaborn as sns
import matplotlib as mpl

def rm_color_cats(c_dict, order, cats):
    if cats:
        keys = c_dict.keys()
        pop_list = []
        for key in keys:
            if key not in cats:
                pop_list.append(key)
        for p in pop_list:
            del c_dict[p]
        order = [o for o in order if o in cats]
    return c_dict, order

def get_talon_nov_colors(cats=None):
    c_dict = {'Known': '#009E73',
              'ISM': '#0072B2',
              'ISM_rescue': '#0072B2',
              'NIC': '#D55E00',
              'NNC': '#E69F00',
              'Antisense': '#000000',
              'Intergenic': '#CC79A7',
              'Genomic': '#F0E442'}
    order = ['Known', 'ISM', 'ISM_rescue', 'NIC', 'NNC', 'Antisense', 'Intergenic', 'Genomic']

    c_dict, order = rm_color_cats(c_dict, order, cats)
    return c_dict, order

def get_pop_colors():
    c_dict = {'AJI': '#46bff0',
            'HAC': '#4cb33e',
            'ITU': '#db72f2',
            'LWK': '#A09136',
            'MPC': '#eb9d0c',
            'PEL': '#ff3a33',
            'YRI': '#DFBD00',
            'CEU': '#347eed'}
    return c_dict, list(c_dict.keys())


def init_plot_settings(font_scale=2,
                       aspect='square',
                       subplot_r=None,
                       subplot_c=None,
                       w=None,
                       h=None):
    """
    Initialize default plotting settings
    """
    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42

    if aspect=='square':
        mpl.rcParams['figure.figsize'] = (5,5)
    elif aspect=='rectangle':
        mpl.rcParams['figure.figsize'] = (7,5)
    
    if w and h:
        mpl.rcParams['figure.figsize'] = (w,h)
    if subplot_r and subplot_c:
        mpl.rcParams['figure.figsize'] = (5*subplot_c, 5*subplot_r)
    