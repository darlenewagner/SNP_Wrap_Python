import random
from datetime import datetime, timedelta

min_year=2121
max_year=2121

start = datetime(min_year, 6, 1, 00, 00, 00)
years = max_year - min_year + 1
end = start + timedelta(days=365*years)

for i in range(17):
    random_date = start + (end - start)*random.random()
    #print datetime object as date only
    print(random_date.date())

#done
