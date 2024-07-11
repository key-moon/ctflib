from ..general import *

def pause():
  input("[+] press any key to continue...")

TOKEN = "eyJhbGciOiJFZERTQSIsInR5cCI6IkpXVCJ9.eyJleHAiOjE3MDAwNDI0MDAsImlhdCI6MTY5OTkyNDgwMCwidGVhbUlkIjo4LCJyb2xlIjoxfQ.LhJ6YSVmcdlebhcnU53oRF2jQ-JRkNi3hyAD0Aodzc1vc6lKwriqjc9FV2ZnL3b8dEY96k7bfjbhSfhy9MU3Aw"
def submit_flag(flag):
    print(requests.post(
        'https://final2023.hitconctf.com/v2/team/me/flag',
        json={ "flag": flag },
        headers={ "Authorization": TOKEN },
    ).text)
