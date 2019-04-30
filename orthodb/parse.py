import json
import requests

req = "https://www.orthodb.org/siblings?id=716at7742&limit=2"

response = requests.get(req)
json_data = json.loads(response.text)


print(json_data["data"][0].keys())