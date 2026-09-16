import requests
import pandas as pd

# 1. 取得したい場所（マサチューセッツ州 Westbrook付近）の緯度・経度
latitude = 42.431868
longitude = -72.664232


# 2. 取得したい期間を設定（例：2023年1月1日〜2023年12月31日）
start_date = "1997-01-01"
end_date = "2025-12-31"

# Open-MeteoのAPIエンドポイントとパラメータ設定
url = "https://archive-api.open-meteo.com/v1/archive"
params = {
    "latitude": latitude,
    "longitude": longitude,
    "start_date": start_date,
    "end_date": end_date,
    "daily": "precipitation_sum",  # 1日の総降水量
    "timezone": "America/New_York" # タイムゾーン
}

print("データを取得中...")
response = requests.get(url, params=params)

# データが正常に取得できたか確認
if response.status_code == 200:
    data = response.json()

    # 3. データをPandasのデータフレーム（表形式）に変換
    df = pd.DataFrame({
        "Date": data["daily"]["time"],
        "Precipitation_mm": data["daily"]["precipitation_sum"]
    })

    # 4. CSVファイルとして保存（ファイル名にMAを追加しています）
    save_path = r"C:\Users\shoah\codes\Westbrook_MA_precipitation.csv"
    df.to_csv(save_path, index=False, encoding="utf-8-sig")

    print(f"成功: CSVファイルの保存が完了しました！\n保存先: {save_path}")
    
    # 最初の5行だけ画面にも表示して確認
    print("\n【取得データの一部】")
    print(df.head())
    
else:
    print(f"エラーが発生しました。ステータスコード: {response.status_code}")