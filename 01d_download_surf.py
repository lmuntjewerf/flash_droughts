import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels-monthly-means',
    {
        'product_type': 'monthly_averaged_reanalysis',
        'variable': [
            'geopotential', 'land_sea_mask', 'soil_type',
        ],
        'year': '2000',
        'month': '01',
        'time': '00:00',
        'area': [
            72, -15, 34,
            41,
        ],
        'format': 'netcdf',
    },
    'surf_para.nc')

