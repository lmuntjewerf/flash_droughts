README

### TESTING to make a flash drought definition

1. Rijnbasin 1D timeserie maken voor precip en PET dagelijks & maandelijks
    Get 2 members h010 s010
    SPEI1 berekenen dagelijks voor 1 dag, bekijken, checken
    SPEI1 montly berekenen
    SPEI1 berekenen dagelijks voor alle dagen?
    - Wat is SPEI1, is dat 30 dagen herhaalpatroon? 
    - Op welke dagen starten? Elke 1e vd maand? 
    deltaSPEI1 maken


Get h010 and s010 pet & precip daily & montly. 
`for letter in h s; do for freq in day Amon; do for var in pet pr; do ecp ec:/nklm/LENTIS/ec-earth/cmorised_by_var/${letter}xxx/${freq}/${var}/${letter}010_${freq}_${var}.tar.gz .; done ; done; done`

2. Rijnbasin 1D timeserie maken voor mrsos and mrso dagelijks & maandelijks
    Get 2 members h010 s010
    - hoe ziet het bodemvocht eruit in/na de dagen met lage SPEI?
    - hoe ziet het bodemvocht eruit na een snelle SPEI decline
    - hoe snel reageert het bodemvocht?

